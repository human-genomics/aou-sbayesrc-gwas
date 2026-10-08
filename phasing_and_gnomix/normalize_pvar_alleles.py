"""Right-trim REF/ALT alleles in a plink2 pvar (after multi-allelic split).

When plink2 splits a multi-allelic VCF/pgen site with `--make-pgen
'multiallelics=-'`, it keeps the original REF as the anchor for every split
row. If one of the original ALT alleles was an indel, the kept REF can be
longer than necessary for the other (SNP) ALT alleles. The result is rows
like ``REF=TATG ALT=CATG`` for what is actually a single T→C SNP. AoU's
plink_bed shipping format and our SBayesRC alignment table both use the
minimal (parsimonious) representation T:C, so plink2's --extract against
``chr:pos:T:C`` misses these rows by ~3% of SBayesRC variants on chr22.

Right-trimming the shared suffix of REF and ALT (keeping at least one base
in each) recovers the minimal form for affected SNPs and for indels that
plink2's split left over-padded on the right. It is purely a textual change
to the pvar — the pgen indexes variants by row order and uses allele indices
(0=REF, 1=ALT), which are unaffected by relabeling. Genotypes do not move.

We do not left-trim, because that would require shifting POS, which in turn
would require touching INFO fields whose values may depend on position. Left
alignment is already enforced by the AoU upstream pipeline.

Usage:
    python3 normalize_pvar_alleles.py --in chr22.split.pvar \
                                       --out chr22.split.pvar.normalized
or in-place via the orchestrator's --in == --out.
"""

import argparse
import sys
from pathlib import Path


def right_trim(ref: str, alt: str) -> tuple[str, str]:
    """Drop the largest common suffix while keeping >=1 base each."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    return ref, alt


def normalize(in_path: Path, out_path: Path) -> dict:
    rewritten = 0
    total = 0
    with in_path.open() as fin, out_path.open("w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                fout.write(line)
                continue
            total += 1
            chrom, pos, vid, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            new_ref, new_alt = right_trim(ref, alt)
            if (new_ref, new_alt) != (ref, alt):
                rewritten += 1
                fields[3] = new_ref
                fields[4] = new_alt
                # Update ID if it follows the chr:pos:REF:ALT convention.
                # We rebuild from the (already-correct) chrom + pos so we don't
                # rely on parsing the old ID.
                fields[2] = f"{chrom}:{pos}:{new_ref}:{new_alt}"
            fout.write("\t".join(fields) + "\n")
    return {"total": total, "rewritten": rewritten}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in", dest="in_path", type=Path, required=True)
    parser.add_argument("--out", dest="out_path", type=Path, required=True)
    args = parser.parse_args()

    if args.in_path == args.out_path:
        # in-place: write to .tmp then rename
        tmp = args.out_path.with_suffix(args.out_path.suffix + ".normalize.tmp")
        stats = normalize(args.in_path, tmp)
        tmp.replace(args.out_path)
    else:
        stats = normalize(args.in_path, args.out_path)

    print(
        f"normalize_pvar_alleles: total={stats['total']} "
        f"rewritten={stats['rewritten']} "
        f"({100.0 * stats['rewritten'] / max(stats['total'], 1):.2f}%)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
