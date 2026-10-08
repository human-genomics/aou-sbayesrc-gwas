"""One full-cohort chromosome: extract, phase, verify, and shard once."""
from __future__ import annotations

from contextlib import ExitStack
from pathlib import Path
import hashlib
import resource
import shutil
import time
import numpy as np
import pandas as pd
from phasing_and_gnomix.common import run, read_psam, write_json, complete, log
from phasing_and_gnomix.extract import prepare, pvar_records


def check_and_shard(qc, phased, samples, output, shard_size):
    import pgenlib
    qc, phased, output = map(Path, (qc, phased, output))
    q_ids, p_ids = read_psam(qc.with_suffix('.psam')), read_psam(phased.with_suffix('.psam'))
    if list(q_ids.IID) != list(samples.IID) or list(p_ids.IID) != list(samples.IID):
        raise ValueError('Phasing changed the cohort or sample order')
    variants = list(pvar_records(qc.with_suffix('.pvar')))
    if variants != list(pvar_records(phased.with_suffix('.pvar'))):
        raise ValueError('Phasing changed the variant order or alleles')
    n, v = len(samples), len(variants)
    # Explicitly retain original FID/IID mapping after IID-only VCF interchange.
    samples[['FID','IID']].rename(columns={'FID':'#FID'}).to_csv(phased.with_suffix('.psam'), sep='\t', index=False)
    shard_root = output / 'shards'
    shard_root.mkdir(parents=True, exist_ok=True)
    n_shards = (n + shard_size - 1) // shard_size
    writers = []
    expected_digests = [hashlib.sha256() for _ in range(n_shards)]
    imputed = observed = 0
    block = 64  # ~150 MB per allele buffer at the production cohort size.
    with ExitStack() as stack:
        src = stack.enter_context(pgenlib.PgenReader(str(qc.with_suffix('.pgen')).encode(), raw_sample_ct=n, variant_ct=v))
        dest = stack.enter_context(pgenlib.PgenReader(str(phased.with_suffix('.pgen')).encode(), raw_sample_ct=n, variant_ct=v))
        for b in range(n_shards):
            s, e = b * shard_size, min((b+1)*shard_size, n)
            prefix = shard_root / f'batch_{b:05d}'
            samples.iloc[s:e][['FID','IID']].rename(columns={'FID':'#FID'}).to_csv(prefix.with_suffix('.psam'), sep='\t', index=False)
            writers.append(stack.enter_context(pgenlib.PgenWriter(str(prefix.with_suffix('.pgen')).encode(), e-s, v, nonref_flags=False, allele_ct_limit=2, hardcall_phase_present=True)))
        for start in range(0, v, block):
            end = min(start+block, v)
            a = np.empty((end-start,2*n), np.int32)
            b = np.empty_like(a)
            ph = np.empty((end-start,n), np.uint8)
            src.read_alleles_range(start,end,a)
            dest.read_alleles_and_phasepresent_range(start,end,b,ph)
            if (b < 0).any() or (b > 1).any() or (((b[:,0::2] != b[:,1::2]) & (ph == 0)).any()):
                raise ValueError('Beagle output contains missing, unphased, or non-biallelic calls')
            called = (a[:,0::2] >= 0) & (a[:,1::2] >= 0)
            if (((a[:,0::2]+a[:,1::2]) != (b[:,0::2]+b[:,1::2])) & called).any():
                raise ValueError('Beagle changed an originally called unordered genotype')
            observed += int(called.sum())
            imputed += int((~called).sum())
            for shard, writer in enumerate(writers):
                left, right = shard*shard_size, min((shard+1)*shard_size,n)
                chunk = np.ascontiguousarray(b[:,2*left:2*right])
                expected_digests[shard].update(chunk.tobytes())
                writer.append_alleles_batch(chunk, all_phased=True)
            if start % (block*1000) == 0:
                log(f'Verified and sharded {end}/{v} variants')
    # Read each shard once against the streaming digest of its source alleles.
    # Never rescan the huge cohort-wide pfile once per sample shard.
    for shard in range(n_shards):
        left, right = shard*shard_size, min((shard+1)*shard_size,n)
        actual = hashlib.sha256()
        with pgenlib.PgenReader(str(shard_root / f'batch_{shard:05d}.pgen').encode(), raw_sample_ct=right-left, variant_ct=v) as dst:
            for start in range(0,v,4096):
                end=min(start+4096,v)
                b=np.empty((end-start,2*(right-left)),np.int32)
                ph=np.empty((end-start,right-left),np.uint8)
                dst.read_alleles_and_phasepresent_range(start,end,b,ph)
                if (((b[:,0::2] != b[:,1::2]) & (ph == 0))).any():
                    raise ValueError('Shard lost phase information')
                actual.update(b.tobytes())
        if actual.digest() != expected_digests[shard].digest():
            raise ValueError('Shard phase differs from retained Beagle phase')
    return {'samples':n,'variants':v,'shards':n_shards,'observed_genotypes_verified':observed,'missing_genotypes_filled_by_beagle':imputed}


def run_phase(args, config, manifest):
    started = time.time()
    rt, out, scratch = map(Path,(args.runtime,args.output,args.scratch))
    scratch.mkdir(parents=True,exist_ok=True)
    out.mkdir(parents=True,exist_ok=True)
    samples = pd.read_csv(args.samples, sep='\t', dtype={'IID':str,'FID':str})
    qc = prepare(args.chrom, args.source, args.psam, args.keep, args.meta, out, scratch, rt/'plink2', config)
    if len(read_psam(qc.with_suffix('.psam'))) != len(samples):
        raise ValueError('Extraction dropped cohort samples')
    if args.chrom == 22:
        from phasing_and_gnomix.validation import demo
        demo(rt,args.model,config,out/'aou_mask_validation.json',out/'match.npz')
    # Source and split files are disposable copies localized on this Batch worker.
    for prefix in (Path(args.source),scratch/'split'):
        for ext in ('.pgen','.pvar','.psam'):
            path=prefix.with_suffix(ext)
            if path.exists() and path.resolve().is_relative_to(Path('/mnt/data')):
                if path.is_symlink():
                    path.resolve().unlink()
                path.unlink()
    plink_args=[rt/'plink2','--seed',str(config['seeds']['plink'])]
    run([*plink_args,'--pfile',qc,'--export','vcf','bgz','id-paste=iid','--output-chr','chrM','--out',scratch/'input'])
    run([rt/'jre/bin/java',f'-Xmx{args.heap_gb}g','-jar',rt/'beagle.jar',f'gt={scratch}/input.vcf.gz',
         f'map={args.map}',f'out={scratch}/beagle','impute=false',f"seed={config['seeds']['beagle']}",
         f"nthreads={config['beagle_threads']}"])
    phased=out/f'chr{args.chrom}'
    run([*plink_args,'--vcf',scratch/'beagle.vcf.gz','--make-pgen','--no-pheno','--output-chr','chrM','--out',phased])
    shutil.copyfile(scratch/'beagle.log',out/'beagle.log')
    stats=check_and_shard(qc,phased,samples,out,config['shard_size'])
    stats.update(wall_seconds=time.time()-started,seed=config['seeds']['beagle'],heap_gb=args.heap_gb,
                 beagle_version=config['beagle_version'],gnofix_required_for_every_sample=True,
                 peak_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    write_json(out/'phasing.json',stats)
    complete(out,manifest['manifest_id'],{'stage':'phase','chrom':args.chrom,**stats})
