# AoU union phasing and Gnomix/Gnofix

Run this optional workflow after `bash get_genotypes.sh` has produced its
ADMIXTURE classification and PCA-fit keep list. It selects everyone outside the
saved European classification plus a seeded random sample of 50,000 PCA-fit
Europeans, phases the **union of pretrained Gnomix and SBayesRC hg38 variants**,
and runs Gnomix and Gnofix for every selected person on all 22 autosomes.

All scripts for this workflow are in `phasing_and_gnomix/`. Public Git contains
code, configuration, seeds and synthetic tests. Private working files default to
the repository's ignored `data/phasing_and_gnomix/`; large genotypes and results
are stored in the workspace's GCS bucket and temporary Batch-worker disks.
No karyograms or other participant plots are generated.

## Inputs from the main pipeline

The workflow reads these files below the main pipeline output prefix, normally
`${WORKSPACE_BUCKET_URI}/sbayesrc_genotypes/`:

| File | Purpose |
|---|---|
| `europeans/classified_european_iids.txt` | Saved European classification; everyone outside it is selected |
| `pca_eur/fit_pca_iids.txt` | Unrelated European candidates for the reproducible 50,000-person sample |
| `statgen/aou_admixture_k6.tsv` | Verify the keep list against the configured K=6 classification rule |

The **genotype source is the original AoU `acaf_threshold/pgen` callset** for the
same data release. The main pipeline's `wgs_pfiles/` and GWAS Step 2 pfiles
already contain only SBayesRC sites, so they cannot supply all Gnomix features.
The coordinator downloads small sample metadata files; it does not download
the full participant genotype files to the notebook.

The default classification check matches `classify_admixture_europeans.py`:
European >=80%, with African, American, East Asian and Oceanian each <=10%.
South Asian has no separate ceiling in this rule. If the main pipeline used
different thresholds, supply a matching private configuration to `prepare`.
PCA-fit IDs must be contained in the saved European IDs, and all lists must
refer to the same WGS sample universe. Preparation fails on inconsistencies.

Public references are downloaded automatically and frozen in a run manifest:

- Original [pretrained Gnomix models](https://github.com/AI-sandbox/gnomix),
  with their original eight ancestry labels and feature order.
- The pinned [gnomix-1000g adapter](https://github.com/human-genomics/gnomix-1000g)
  and UCSC hg19-to-hg38 chain for model-coordinate/allele liftover.
- [`sbayesrc_hg38.csv`, release v1.0](https://github.com/human-genomics/sbayesrc-liftover/releases/tag/v1.0).
  This table is already hg38 and is not lifted again.
- Pinned GRCh38 maps from GLIMPSE, converted to Beagle's four-column PLINK map
  format, Beagle 5.5 build `27Feb25.75f`, Java 17, and PLINK 2 build 29 Jan 2025.

The model's original alleles are audited before any one-character conversion.
This pinned model is SNV-only. The statistical model is not retrained: hg38
genotypes are mapped back to its original feature order and allele coding.

## Run from an AoU Jupyter terminal

Use an authorized Controlled Tier Workbench Jupyter environment with access to
the original source data, a writable workspace bucket, Google Batch, and the
workspace pet service account. The notebook control environment requires
Python >=3.10. Model inference uses a separate pinned Python 3.9 runtime.

```bash
git clone https://github.com/human-genomics/aou-sbayesrc-gwas.git
cd aou-sbayesrc-gwas
python -m pip install -r phasing_and_gnomix/requirements-orchestrator.txt

# After get_genotypes.sh has produced the classification and PCA-fit inputs:
bash phasing_and_gnomix/run.sh prepare --run-id union-v1
bash phasing_and_gnomix/run.sh dry-run --run-id union-v1

# Starts paid Batch work. The notebook coordinator runs detached.
nohup bash phasing_and_gnomix/run.sh run --run-id union-v1 \
  > data/phasing_and_gnomix/runs/union-v1/coordinator.log 2>&1 &

bash phasing_and_gnomix/run.sh status --run-id union-v1
```

For an existing checkout, start with dependency installation; do not clone it
again. `prepare` downloads tools/references, runs public/synthetic validation,
checks the cohort and source objects, and stages immutable inputs. It submits
no Batch jobs. `dry-run` prints the chr22 worker request without submitting it.
`run` first executes a public/synthetic cloud smoke test, then the entire chr22
cohort through phasing, Gnomix/Gnofix, corrected scoring export and tract audit.
Only after every chr22 shard passes does it submit the remaining chromosomes.
There is no additional phasing concurrency cap; Google quota and VM capacity
may queue requests. Inference begins as each chromosome finishes phasing.

Closing the browser or terminal does not stop a `nohup` coordinator. **Keep
the Jupyter environment running**: already submitted Batch jobs are independent,
but subsequent inference submission and final aggregation need the coordinator.
To resume after it stops, repeat the same `run` command with the same `--work`
and `--run-id`. An exclusive lock prevents two coordinators from owning a run.

The workspace project/account and bucket are discovered from the notebook,
preferring the mounted bucket over a potentially stale `WORKSPACE_BUCKET`.
Paths and compute placement can be specified explicitly:

```bash
bash phasing_and_gnomix/run.sh prepare --run-id union-v1 \
  --project YOUR_WORKSPACE_PROJECT \
  --workspace-bucket gs://YOUR_WORKSPACE_BUCKET \
  --pipeline-uri gs://YOUR_WORKSPACE_BUCKET/sbayesrc_genotypes \
  --data-version v9 \
  --resources data/my_union_resources.json
```

Other overrides include `--source-uri`, `--output-uri`, `--region`,
`--service-account`, `--network`, `--subnetwork`, `--config`, `--seeds` and
`--work`. Specify input/resource overrides during preparation; run/resume uses
the frozen manifest. Use a new run ID/output prefix when changing inputs,
settings or code. Outputs from different variant sets cannot be mixed.

## Compute and reproducibility

`resources.json` contains explicit RAM, Java heap, thread count and machine
type for every chromosome. It records the current full-cohort v9 allocations,
including large M1/M3 machines. **Validation of the larger union chromosomes
is still in progress.** These settings are a starting profile for a similar
cohort, not a universal memory guarantee. A successful chr22 pilot verifies
the end-to-end workflow but does not establish memory requirements for larger
or denser chromosomes. Review the resource file before starting another cohort.

The scientific Beagle parameters remain at their defaults. The command sets
input/output/map paths, `impute=false`, the tracked seed, thread count and Java
heap. It uses no external reference panel or pedigree. Missing genotypes at
retained sites are filled by Beagle as part of population phasing. Phasing
always includes the whole selected cohort; later sample sharding is solely
for independent Gnomix/Gnofix inference and output storage.

`seeds.json` records cohort, PLINK, Beagle, inference and validation seeds.
The cohort RNG uses pinned NumPy/PCG64 over sorted eligible IDs. Each Gnofix
sample seed is derived from the recorded seed, chromosome and IID with SHA-256,
so inference order and shard boundaries do not change it. The run also freezes
source object generations, public checksums, tools, model audit, resource
settings, sample order and source code. Bitwise identity across different
processors and numerical runtimes is not guaranteed.

Temporary API submission failures are reconciled against the exact attempted
Batch ID and run/unit labels. Accepted jobs are adopted on resume. Retrying an
unaccepted request requires an authoritative 404, a grace period and a second
absence check. An uncertain lookup never authorizes a replacement. Failed paid
compute is reported and is not automatically rerun with new memory or scientific
settings; unrelated chromosomes can continue. Worker stdout/stderr is in GCS.

## Variant handling

The production extraction follows the established seven-command path: candidate
selection, multiallelic split, missingness measurement, sorted union selection,
one VCF export, Beagle, and phased-PGEN import. It does not assemble separate
model/scoring genotype copies or concatenate/sort a second VCF.

1. Select source records at lifted model or SBayesRC positions, split
   multiallelics, and assign `chr:pos:REF:ALT` IDs.
2. Right-trim padded split alleles using the same logic as the repository's
   `normalize_pvar_alleles.py`, before applying the SNV restriction.
3. Match position **and alleles**, accounting for reference/alternate reversal
   and model liftover strand. Reject ambiguous source/model mappings.
4. Remove variants with missingness **above 10%** across the selected cohort;
   exactly 10% passes. No extra MAF, HWE, FILTER/QUAL/SCORE/
   CALIBRATION_SENSITIVITY or sample `--mind` filter is applied. The source
   already has AoU's ACAF threshold selection.
5. Phase the retained union and verify observed unordered genotypes, full
   phase, sample/variant order, and shard readback against the source.

Gnomix consumes only its expected model features. Additional SBayesRC variants
do not enter the classifier. Absent, failed-liftover or QC-excluded model
features are encoded as model-reference (`0`) for inference, following the
pinned adapter and public-demo validation. This is a modeling convention, not
a claim that an absent participant genotype is biologically homozygous
reference. Unavailable scoring SNPs are not fabricated in the output PGEN.

`feature_qc.tsv.gz`, `sbayesrc_feature_qc.tsv.gz` and
`variant_missingness.tsv.gz` record matching and QC outcomes. A retained-model
coverage floor of 80% stops a chromosome with insufficient usable features.

## GCS output layout

Default root: `${PIPELINE_URI}/phasing_and_gnomix/${RUN_ID}/`.

```text
inputs/                         frozen cohort, code, runtime, models, maps, scoring lists
preflight/cloud_smoke/          public/synthetic worker validation
phased/chrN/
  chrN.{pgen,pvar,psam}          full retained UNION, with Beagle phase
  beagle.log, phasing.json       phasing settings, runtime and peak child RSS
  match.npz, windows.tsv         trained model-feature and window mappings
  scoring_variants.tsv.gz        scoring index -> source variant, alleles, rsid, window
  *_qc.tsv.gz                   matching/exclusion audits
  variant_missingness.tsv.gz
  shards/batch_XXXXX.{pgen,psam}  union genotypes split by samples for inference
  COMPLETE.json                 inventory, hashes and verification results
inference/chrN/batch_XXXXX/
  scoring_genotypes_gnofix.{pgen,pvar,psam}
  scoring_manifest.json
  tracts_gnofix.tsv.gz
  parts/part_XXXXX.npz
  gnofix_switches.tsv.gz, sample_status.tsv.gz
  ancestry_totals.npz, tract_audit.json
  COMPLETE.json
results/
  global_ancestry_gnofix.tsv.gz
  sample_index.tsv.gz
  chromosome_index.tsv.gz
  scoring_file_index.tsv.gz
  tract_file_index.tsv.gz
  DATASET_VERSION.json, COMPLETE.json
logs/dsub/                      per-job stdout/stderr
state.json                      mutable submission/completion progress
```

`scoring_genotypes_gnofix` contains **retained SBayesRC SNPs with actual Gnofix
phase correction applied**, stored by chromosome and sample shard. The full
union `phased/chrN/chrN.pgen` has Beagle-only phase. Both are labeled explicitly;
do not substitute the latter for corrected scoring genotypes. Saved swap
matrices also describe the Gnofix corrections on the model windows.

`scoring_variants.tsv.gz` gives each scoring PGEN row's `scoring_variant_index`,
union `source_variant_index`, rsid, source REF/ALT, SBayesRC allele reversal and
zero-based model `window`. A matched model site uses its exact model window;
an additional SBayesRC site uses the nearest reliably lifted model feature,
with the upstream feature winning a distance tie. That same window determines
the genotype swap, ancestry and posterior. Additional sites do not receive
independent SNP-level ancestry estimates.

For a sample in `parts/part_XXXXX.npz`, find its index `i` in `samples`.
`labels[2*i, window]` and `labels[2*i+1, window]` describe corrected haplotypes
A and B; `called_posterior` supplies the probability of each called ancestry.
These files store called-ancestry support, not the full eight-class probability
vector. Label codes follow `config.json`'s ancestry order. Haplotype A/B does
not imply maternal/paternal origin. Align predictor effect alleles to the source
alleles recorded in the scoring map before scoring.

Tracts have one row per consecutive same-ancestry run on a haplotype:
`sample, haplotype, chrom, start_hg38, end_hg38, start_hg19, end_hg19,
start_cM, end_cM, ancestry, n_windows, mean_posterior`.
Unreliable lifted tract endpoints use `-1`. All eight labels remain separate:
AFR, AHG, EUR, WAS, SAS, EAS, NAT and OCE. EUR/WAS and AFR/AHG are not merged.
Genome-wide percentages are weighted by genetic window length across both
haplotypes and all 22 autosomes. Final aggregation requires every expected
sample/chromosome/shard and verified scoring corrections.

## Validation and publication

The public demo checks agreement with the pinned 1000G adapter, batch-invariant
Gnofix calls/swaps, sensitivity to absent model features, tract/posterior
consistency, and real corrected scoring PGEN readback. Synthetic integration
tests run actual PLINK and Beagle on a union containing padded multiallelics,
allele reversals, duplicate sites, extra scoring variants and missing calls.
Coordinator tests cover pilot gating, API reconciliation, immutable input
boundaries and rejection of incomplete/mixed results.

```bash
# After prepare has built the pinned runtime:
export LD_LIBRARY_PATH="$PWD/data/phasing_and_gnomix/runtime/libgomp/lib:$PWD/data/phasing_and_gnomix/runtime/libstdcxx/lib:$PWD/data/phasing_and_gnomix/runtime/libgcc/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
data/phasing_and_gnomix/runtime/python/bin/python3 -m pytest -q phasing_and_gnomix/tests
# The fake-cloud coordinator exercise uses the notebook's Google client libraries:
python -m phasing_and_gnomix.tests.exercise_scheduler
```

The portable entry point is tested locally with synthetic/public data. A fresh
all-22 AoU launch through this new entry point has not yet been completed; the
ongoing cohort run uses the preceding frozen worker packages. The reusable
scientific modules were derived from that union workflow, without changing its
extraction/QC, Beagle defaults, model inputs or per-person Gnofix algorithm.
Check resource requirements and workspace quotas for each new environment.

Keep all generated manifests, sample lists, logs, genotypes, ancestry outputs
and run-specific notes in ignored private storage. The source-only worker bundle
uses an explicit file allowlist. No release upload, Git commit or push is
performed by this workflow.
