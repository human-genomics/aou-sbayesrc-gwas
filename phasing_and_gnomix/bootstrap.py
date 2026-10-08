"""Build a relocatable, checksummed runtime and fetch public references only."""
from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import urllib.request
import zipfile

from .common import read_json, run, sha256, write_json, log

PUBLIC = {
    'plink2.zip': ('https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250129.zip', '5477728b725f2aedd6fab33139dc884c215dd8dcb58583cadf373d3dc50ba8b7'),
    'adapter.tar.gz': ('https://github.com/human-genomics/gnomix-1000g/archive/0472e5b019f359f9d5b2084e06b0a85973ecba6a.tar.gz', 'ee29b46049154f88eba0f1bb348b775fb66448b7ca2f53e3a13ffb5bc4831835'),
    'python.tar.gz': ('https://github.com/astral-sh/python-build-standalone/releases/download/20241016/cpython-3.9.20%2B20241016-x86_64-unknown-linux-gnu-install_only.tar.gz', 'c20ee831f7f46c58fa57919b75a40eb2b6a31e03fd29aaa4e8dab4b9c4b60d5d'),
    'jre.tar.gz': ('https://github.com/adoptium/temurin17-binaries/releases/download/jdk-17.0.20.1%2B1/OpenJDK17U-jre_x64_linux_hotspot_17.0.20.1_1.tar.gz', '0b2b640e3046b64c8ec504de0ab9d91bb5610182bda21fad454681ce54d45a62'),
    'beagle.jar': ('https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar', '7319f4af9638be05c18dcc1bfb8fb41a58a09293507ebf0d54617d0e40df5a70'),
    'gnomix.tar.gz': ('https://github.com/AI-sandbox/gnomix/archive/cd15f65eadfca8e2bbf209eb69a82978a200c5fd.tar.gz', '846ef1ba78b569c04acaadf1c6c6e6698de881b005461f9d3ea70bf29bd6b6d8'),
    'pretrained_gnomix_models.tar.gz': ('https://drive.usercontent.google.com/download?id=1Q0zg9zqTaZUvt42uxzE_0gcfnFvjwi33&export=download&confirm=t', '06cc751d06ff1244dce7ffe6366e27cbedd49992f79d37b346faeed4ea748ff1'),
    'hg19ToHg38.over.chain.gz': ('https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz', '5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d'),
    'libgomp.tar.bz2': ('https://conda.anaconda.org/conda-forge/linux-64/libgomp-12.2.0-h65d4601_19.tar.bz2', '81a76d20cfdee9fe0728b93ef057ba93494fd1450d42bc3717af4e468235661e'),
    'libstdcxx.tar.bz2': ('https://conda.anaconda.org/conda-forge/linux-64/libstdcxx-ng-12.2.0-h46fd767_19.tar.bz2', '0289e6a7b9a5249161a3967909e12dcfb4ab4475cdede984635d3fb65c606f08'),
    'libgcc.tar.bz2': ('https://conda.anaconda.org/conda-forge/linux-64/libgcc-ng-12.2.0-h65d4601_19.tar.bz2', 'f3899c26824cee023f1e360bd0859b0e149e2b3e8b1668bc6dd04bfc70dcd659')
}


def fetch(url, path, expected=None):
    path = Path(path)
    # A path alone cannot identify an unchecksummed URL (e.g. a map commit
    # changed in the config). Whole-reference cache reuse is checked by prepare.
    if path.exists() and expected is not None and sha256(path) == expected:
        return sha256(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    log('Downloading ' + path.name)
    tmp = path.with_suffix(path.suffix + '.part')
    try:
        run(['curl', '-fL', '--connect-timeout', '20', '--retry', '4', '--retry-delay', '3', '--silent', '--show-error', url, '-o', tmp])
    except subprocess.CalledProcessError:
        if path.name != 'hg19ToHg38.over.chain.gz':
            raise
        # Byte-identical mirror; the original UCSC checksum remains mandatory.
        mirror = 'https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/over.chain/hg19ToHg38.over.chain.gz'
        run(['curl','-fL','--retry','3','--silent','--show-error',mirror,'-o',tmp])
    digest = sha256(tmp)
    if expected and digest != expected:
        raise ValueError('Checksum mismatch for ' + path.name)
    tmp.replace(path)
    return digest


def unpack(archive, dest, strip=False):
    dest = Path(dest)
    if (dest / '.unpacked').exists():
        if (dest / '.unpacked').read_text().strip() != sha256(archive):
            raise ValueError('Cached extraction differs from the pinned archive; use a fresh work directory')
        return
    dest.mkdir(parents=True, exist_ok=True)
    run(['tar', '-xf', archive, '-C', dest, *(['--strip-components=1'] if strip else [])])
    (dest / '.unpacked').write_text(sha256(archive) + '\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--config', type=Path, default=Path(__file__).with_name('config.json'))
    args = p.parse_args()
    cfg = read_json(args.config)
    if cfg['gnomix_commit'] not in PUBLIC['gnomix.tar.gz'][0] or cfg['beagle_version'] not in PUBLIC['beagle.jar'][0]:
        raise ValueError('Configured software differs from the pinned archives')
    dl, rt = args.work / 'downloads', args.work / 'runtime'
    dl.mkdir(parents=True, exist_ok=True)
    rt.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(4) as ex:
        futures = [ex.submit(fetch, u, dl / n, h) for n, (u, h) in PUBLIC.items()]
        for f in futures:
            f.result()
    if cfg['adapter_commit'] not in PUBLIC['adapter.tar.gz'][0]:
        raise ValueError('Adapter commit differs from the pinned archive')
    fetch(cfg['sbayesrc_url'], dl / 'sbayesrc_hg38.csv', cfg['sbayesrc_sha256'])
    unpack(dl / 'python.tar.gz', rt)
    unpack(dl / 'jre.tar.gz', rt / 'jre', True)
    unpack(dl / 'gnomix.tar.gz', rt / 'gnomix', True)
    unpack(dl / 'adapter.tar.gz', rt / 'adapter', True)
    unpack(dl / 'libgomp.tar.bz2', rt / 'libgomp')
    unpack(dl / 'libstdcxx.tar.bz2', rt / 'libstdcxx')
    unpack(dl / 'libgcc.tar.bz2', rt / 'libgcc')
    # Conda normally creates the SONAME link via _openmp_mutex. This runtime
    # extracts packages directly, so provide it explicitly for minimal images.
    soname=rt/'libgomp/lib/libgomp.so.1'
    if not soname.exists():soname.symlink_to('libgomp.so.1.0.0')
    shutil.copyfile(dl / 'beagle.jar', rt / 'beagle.jar')
    with zipfile.ZipFile(dl / 'plink2.zip') as archive:
        (rt / 'plink2').write_bytes(archive.read('plink2'))
    if sha256(rt / 'plink2') != '28e45a65689bacf2b6675898abc6fc2eb2a5de19ae21a60ccec0dd37e6b21ce2':
        raise ValueError('Pinned PLINK binary checksum differs')
    (rt / 'plink2').chmod(0o755)
    req = Path(__file__).with_name('requirements.txt')
    if not (rt / 'requirements.sha256').exists() or (rt / 'requirements.sha256').read_text().strip() != sha256(req):
        run([rt / 'python/bin/python3', '-m', 'pip', 'install', '--disable-pip-version-check', '-r', req, 'pytest==8.3.5'])
        run([rt / 'python/bin/python3', '-m', 'pip', 'freeze'], stdout=open(rt / 'requirements.lock.txt', 'w'))
        (rt / 'requirements.sha256').write_text(sha256(req) + '\n')
    maps = args.work / 'references/maps'
    maps.mkdir(parents=True, exist_ok=True)
    for c in cfg['chromosomes']:
        url = f"https://raw.githubusercontent.com/odelaneau/GLIMPSE/{cfg['map_commit']}/maps/genetic_maps.b38/chr{c}.b38.gmap.gz"
        fetch(url, maps / f'chr{c}.gmap.gz')
        prev_pos, prev_cm = -1, -1.0
        with gzip.open(maps / f'chr{c}.gmap.gz', 'rt') as src, open(maps / f'chr{c}.map', 'w') as dst:
            assert src.readline().split() == ['pos', 'chr', 'cM']
            for line in src:
                pos, chrom, cm = line.split()
                assert int(chrom) == c and int(pos) > prev_pos and float(cm) >= prev_cm
                dst.write(f'chr{c}\tchr{c}:{pos}\t{cm}\t{pos}\n')
                prev_pos, prev_cm = int(pos), float(cm)
    manifest = {str(p.relative_to(args.work)): sha256(p) for p in sorted(dl.iterdir()) if p.is_file() and not p.name.endswith('.part')}
    manifest.update({str(p.relative_to(args.work)): sha256(p) for p in sorted(maps.iterdir())})
    manifest['runtime/plink2'] = sha256(rt / 'plink2')
    manifest['runtime/requirements.lock.txt'] = sha256(rt / 'requirements.lock.txt')
    write_json(args.work / 'references/public_inputs.json', manifest)
    # A previous runtime archive cannot represent newly installed native libraries.
    (args.work/'runtime.tar.gz').unlink(missing_ok=True)
    log('Public runtime and references ready')


if __name__ == '__main__':
    main()
