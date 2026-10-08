#!/bin/bash
set -euo pipefail
printf '%s  %s\n' "$RUNTIME_SHA256" "$RUNTIME" "$CODE_SHA256" "$CODE" | sha256sum -c -
mkdir -p /mnt/data/gnomix_runtime /mnt/data/gnomix_code "$OUTDIR"
tar -xzf "$RUNTIME" -C /mnt/data/gnomix_runtime
tar -xzf "$CODE" -C /mnt/data/gnomix_code
export PYTHONPATH=/mnt/data/gnomix_code
export LD_LIBRARY_PATH="/mnt/data/gnomix_runtime/libgomp/lib:/mnt/data/gnomix_runtime/libstdcxx/lib:/mnt/data/gnomix_runtime/libgcc/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/gnomix-matplotlib
/mnt/data/gnomix_runtime/jre/bin/java -version
/mnt/data/gnomix_runtime/plink2 --version
/mnt/data/gnomix_runtime/python/bin/python3 -m phasing_and_gnomix.smoke \
  --runtime /mnt/data/gnomix_runtime --model "$MODEL" --meta "$META" \
  --config "$CONFIG" --output "$OUTDIR/smoke_validation.json"
/mnt/data/gnomix_runtime/python/bin/python3 -c 'import os; from phasing_and_gnomix.common import complete,read_json; complete(os.environ["OUTDIR"],read_json(os.environ["MANIFEST"])["manifest_id"],{"stage":"smoke","public_end_to_end_inference_passed":True})'
