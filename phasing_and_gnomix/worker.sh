#!/bin/bash
set -euo pipefail
printf '%s  %s\n' "$RUNTIME_SHA256" "$RUNTIME" "$CODE_SHA256" "$CODE" | sha256sum -c -
mkdir -p /mnt/data/gnomix_runtime /mnt/data/gnomix_code /mnt/data/scratch
tar -xzf "$RUNTIME" -C /mnt/data/gnomix_runtime
tar -xzf "$CODE" -C /mnt/data/gnomix_code
export PYTHONPATH=/mnt/data/gnomix_code
export LD_LIBRARY_PATH="/mnt/data/gnomix_runtime/libgomp/lib:/mnt/data/gnomix_runtime/libstdcxx/lib:/mnt/data/gnomix_runtime/libgcc/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/gnomix-matplotlib
export GNOMIX1000G_TOOLS=/mnt/data/gnomix_runtime
chmod +x /mnt/data/gnomix_runtime/plink2
common=(--runtime /mnt/data/gnomix_runtime --config "$CONFIG" --manifest "$MANIFEST" --output "$OUTDIR" --samples "$SAMPLES" --chrom "$CHROM")
if [[ "$STAGE" == phase ]]; then
    # dsub localizes files by basename; provide a stable prefix explicitly.
    mkdir -p /mnt/data/source
    ln -s "$PGEN" /mnt/data/source/input.pgen
    ln -s "$PVAR" /mnt/data/source/input.pvar
    /mnt/data/gnomix_runtime/python/bin/python3 -m phasing_and_gnomix.worker phase "${common[@]}" \
        --source /mnt/data/source/input --psam "$PSAM" --keep "$KEEP" --meta "$META" \
        --map "$MAP" --model "$MODEL" --heap-gb "$HEAP" --scratch /mnt/data/scratch
else
    /mnt/data/gnomix_runtime/python/bin/python3 -m phasing_and_gnomix.worker infer "${common[@]}" \
        --pgen "$PGEN" --psam "$PSAM" --match "$MATCH" --windows "$WINDOWS" --model "$MODEL" --shard "$SHARD"
fi
