# AoU SBayesRC GWAS pipeline

A reproducible, idempotent pipeline that builds per-chromosome PLINK2 pfiles
holding the ~7.35 M SBayesRC SNPs from the selected All of Us ACAF WGS callset
(`v9` by default), with the variant `ID` column populated with SBayesRC rsids.
It also builds the
~501k direct-SNP bfile, plus a higher-quality direct-SNP bfile filtered with
AoU EUR frequency/missingness metrics. It runs ADMIXTURE K=6 ancestry
projection, writes a European ancestry keep-list, compares against
AoU-provided ancestry calls/fractions when the selected release includes those
inputs, runs KING kinship from the high-quality direct bfile, compares those
estimates to AoU's provided relatedness table, classifies close relationships,
and selects the unrelated European sample set used to fit PCA.
It then builds the PCA-ready SNP bfile from that sample set, fits 20 European
ancestry PCs, and projects those PCs onto all samples in the high-quality
direct bfile. Before GWAS setup, it builds a conservative sample-QC exclusion
list for identical-genotype components of size three or larger and then builds
final REGENIE Step 1/Step 2 genotype inputs with GWAS-specific variant filters.
The final optional example builds a height phenotype/covariate set and runs a
continuous-trait REGENIE GWAS.

Two optional downstream workflows reuse these inputs:
[`pan_aou_gwas/`](pan_aou_gwas/README.md) runs phenotype-wide GWAS, and
[`phasing_and_gnomix/`](phasing_and_gnomix/README.md) prepares phased genotypes
and local ancestry for ancestry-specific scoring. Their scope and prerequisites
are described under [Downstream workflows](#downstream-workflows).

The pipeline **must be run from inside an AoU Verily Jupyter session
terminal** (the standard interactive analysis environment on the All of Us
Researcher Workbench). It cannot be run from an off-platform `wb` CLI
session — VPC SC and IAM together require the orchestrator to run as the
pod's pet service account, which is what the Jupyter pod gives you for free.

From there, the per-chromosome compute runs as **parallel Google Batch jobs
via dsub**, one worker per autosome in `us-central1` (the same region as the
workspace bucket, so output delocalization is intra-region and effectively
line-rate). On a 4-vCPU Jupyter pod a serial run takes ~2–3 days; the dsub
fan-out finishes in ~2 hours.

The pod we developed and tested this on (the minimum-cost default that
worked):

| Component | Value |
|---|---|
| Pod machine type | `n1-highmem-4` |
| vCPUs | 4 |
| RAM | 25 GB |
| Zone | `us-east4-a` (workspace's `terra-default-location`) |
| Local overlay disk | ~492 GB |
| Workspace bucket region | `us-central1` |

The Jupyter pod only orchestrates and submits the dsub tasks — the actual
plink2 work happens on Batch worker VMs we provision in `us-central1`, each
with `--min-cores 4 --min-ram 32 --disk-size 300` by default. You can leave
the Jupyter session open or close the browser (background the run with
`nohup`); the Batch jobs run independently.

By default `get_genotypes.sh` runs against AoU Controlled Tier v9
(`C2025Q4R6`) and writes durable outputs under
`${WORKSPACE_BUCKET}/sbayesrc_genotypes/`, the same output layout used by the
original pipeline. To reproduce the older v8 validation run, set
`AOU_DATA_VERSION=v8` in a clean workspace.

## What this pipeline produces

```
${WORKSPACE_BUCKET}/sbayesrc_genotypes/wgs_pfiles/chr{1..22}.{pgen,pvar,psam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/wgs_pfiles/chr{1..22}.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_pfiles/chr{1..22}.{pgen,pvar,psam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile/chr1_22_merged.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq.filter_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq.sample_missingness_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/aou_admixture_k6.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/europeans/classified_european_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/aou_vs_ours/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/kinship_snp_subset_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/kinship_snp_missingness_threshold_counts.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/aou_hq_direct_rel.kin0
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/qc/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/close_relations.csv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/fit_pca_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/select_pca_europeans.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_ready.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_snp_qc.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_snp_qc.filter_steps.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/aou_projected.sscore
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/aou_pcs.{eigenval,eigenvec,eigenvec.allele}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_eur_counts.acount
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/fit_project_pca.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/genetic_sex/sex_covar.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/genetic_sex/genetic_sex_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/identical_components.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/exclude_identical_component_size_ge3_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/identical_component_sample_qc.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/gwas_genotype_qc.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/gwas_step1_direct.filter_steps.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs_pfiles/chr{1..22}.{pgen,pvar,psam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/gwas_step2_wgs.filter_steps.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/fit_pca_af/gwas_step2_fit_pca_alt_freqs_passing.tsv.gz
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/{phen.txt,covar.txt,training_iids.txt}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/height_gwas.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/step1/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/step2/chr{1..22}/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/lightweight/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/regenie_gwas.summary.tsv

```

Each `summary.tsv` reports `requested / src_variants / src_samples /
biallelic_total / extracted / out_samples / missing / remapped / unmapped`
for that chromosome. A combined summary across all processed chromosomes is
also written to `logs/sbayesrc_extract_summary.tsv` at the end of each run.

## Reporting suppression

README accounting tables suppress exact AoU participant/sample/person-level
counts in the range 1 through 20. Exact zero counts may be reported. This
suppression rule does not apply to thresholds, chromosome numbers, PC counts,
variant/SNP counts, runtime/resource sizing, or other non-participant values.

## Pipeline steps

### Step 1 — `store_sbayesrc_ids.py` (local)

Reads `data/support/sbayesrc_hg38.csv` (chrom, pos, ref, alt, rsid; downloaded
on first run from a public GitHub release) and writes per-chromosome
`data/sbayesrc_ids/chr{N}.extract.txt` (one `chr{N}:pos:REF:ALT` ID per line)
and `chr{N}.idmap.txt` (`chr{N}:pos:REF:ALT \t rsid`). The ID format matches
the AoU `acaf_threshold/plink_bed/chr{N}.bim` second column exactly, so
`plink2 --extract` is a direct ID match.

### Step 2 — `wgs_extract_variants.sh` (dsub fan-out)

Submits **one dsub task per autosome** to Google Batch in `us-central1`.
Each task runs `dsub_extract_worker.sh` on its own Batch worker, with the
per-chromosome AoU pgen/pvar staged in via `dsub --input` and the final
extracted pgen/pvar/psam/summary.tsv delocalized via `dsub --output-recursive`
to `${WORKSPACE_BUCKET}/sbayesrc_genotypes/wgs_pfiles/`. For cdrv9, several
autosomes do not publish their own `.psam` object in `acaf_threshold/pgen`;
the orchestrator verifies that the published autosome `.psam` files are
byte-identical and uses one as a shared sample-file fallback for those
chromosomes.

Per-worker pipeline (four passes, same as the previous serial-on-Jupyter
implementation):

1. **Pass 1** — `plink2 --pfile acaf_threshold.chr{N} --make-pgen 'multiallelics=-' --output-chr chrM --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 10000 --out chr{N}.split`. Splits multi-allelics and assigns `chr:pos:REF:ALT` IDs.

2. **Pass 1.5** (inline awk in `dsub_extract_worker.sh`) — right-trims shared REF/ALT suffix in the split-pvar. plink2's split keeps the longest indel's REF as the anchor for every split row, so SNPs at a multi-allelic site come out over-padded (e.g. `REF=TATG ALT=CATG` for what is really `T:C`). Trimming recovers the minimal form and re-stamps the ID. Purely a text rewrite — pgen indexes by row order with allele codes (0=REF, 1=ALT), so genotypes don't move. On chr22 this raises recovery from 95,130 → 97,960 of 98,065 SBayesRC variants. The canonical, testable reference for this logic is `normalize_pvar_alleles.py`; the worker uses an inline awk version because the default `marketplace.gcr.io/google/ubuntu2204` Batch worker image ships without python3 and VPC SC blocks `apt install`.

3. **Pass 2** — `plink2 --pfile chr{N}.split --extract chr{N}.extract.txt --no-pheno --output-chr chrM --make-pgen --out chr{N}`. Default `pvar-cols` preserves the INFO column (AoU acaf_threshold INFO includes `AC`, `AF`, `AN`, `AS_QUALapprox`, `CALIBRATION_SENSITIVITY`, `QUALapprox`, `SCORE`).

4. **Pass 3** (inline awk) — rewrites the ID column of `chr{N}.pvar` from `chr:pos:REF:ALT` to the SBayesRC rsid via `chr{N}.idmap.txt`.

No QC filters applied — all variants are kept regardless of FILTER/QUAL/
SCORE. ACAF threshold already enforces MAC ≥100 per ancestry. Downstream
code can apply additional filters (e.g. `--extract-if-info
"CALIBRATION_SENSITIVITY < 0.99"`) against the preserved INFO column.

Per-task stdout/stderr lives at
`${WORKSPACE_BUCKET}/sbayesrc_genotypes/logs/dsub/<job-id>.<task>-{stdout,stderr}.log`.

Observed Step 2 accounting from the default cdrv9/v9 clean-workspace run
launched on June 30, 2026:

| Metric | Value |
|---|---:|
| Autosomes processed | 22 |
| Samples per chromosome | 535,662 |
| Requested SBayesRC variants | 7,354,747 |
| Source variants across chromosome pvars | 60,193,204 |
| Biallelic rows after splitting | 133,489,424 |
| Extracted SBayesRC variants | 7,350,153 |
| Missing requested variants | 4,594 |
| Remapped variant IDs | 7,350,153 |
| Unmapped variant IDs | 0 |

### Step 3 — Direct-SNP bfile for REGENIE step 1

`get_genotypes.sh` downloads the UKBB-derived direct-SNP rsid list to
`data/support/direct_snps/ukbb_500k_qc_pass_direct_snps.txt` if it is not
already cached. `prepare_direct_snps.py` cross-references that list against
`sbayesrc_hg38.csv` and the extracted WGS pvars, writing per-chromosome rsid
extract lists and a missing-SNP report for direct SNPs absent from the
AoU/SBayesRC pfiles.

`extract_direct_snps.sh` then builds
`${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_pfiles/chr{N}.{pgen,pvar,psam}`.
Present SNPs are extracted from the WGS pfiles. Direct SNPs absent from the
AoU/SBayesRC pfiles are tracked in `data/direct_snps/missing_direct_snps.tsv`.

`make_direct_bfile.sh` submits one Google Batch/dsub worker in `us-central1`
with an SSD data disk (default: 8 vCPU, 32 GB RAM, 300 GB `pd-ssd`) to merge
the 22 direct pfiles into
`${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile/chr1_22_merged.{bed,bim,fam}`.
The merged bfile contains the direct SNPs present in the AoU extracted pfiles;
absent SNPs are not encoded as all-missing variants. The merged bfile is
intended for REGENIE step 1.

Observed Step 3 accounting from the default cdrv9/v9 run:

| Metric | Value |
|---|---:|
| Requested direct SNPs | 501,477 |
| Present in AoU/SBayesRC WGS pfiles | 501,452 |
| Absent direct SNPs | 25 |
| Merged direct-bfile variants | 501,452 |
| Merged direct-bfile samples | 535,662 |

### Step 4 — Higher-quality direct-SNP bfile

`get_genotypes.sh` downloads `sbayesrc_liftover_results.csv` from the public
SBayesRC liftover release if it is not already cached locally. This file
contains `A1_hg38`, `A2_hg38`, and `A1Freq`, which Step 4 converts to the
SBayesRC hg38 ALT-allele frequency after matching against the AoU direct bfile
REF/ALT alleles.

`make_hq_direct_bfile.sh` builds an AoU EUR keep-list from AoU's computed
ancestry predictions (`ancestry_pred == eur`) and the direct-bfile `.fam`.
It then submits a dsub metrics job that runs plink2 on the raw direct bfile
with that EUR keep-list:

```text
plink2 --bfile chr1_22_merged --keep aou_eur.keep --freq --missing variant-only
```

The local filter builder (`filter_hq_direct_snps.py`) then applies these
filters in order:

```text
1. original UKBB direct-SNP rsid list
2. present in the AoU direct bfile
3. SBayesRC liftover ALT frequency is available and alleles match
4. abs(AoU EUR ALT frequency - SBayesRC ALT frequency) <= 0.04
5. AoU EUR MAF >= 0.007
6. AoU EUR variant missingness <= 0.05
```

The ordered count log is written to:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq.filter_summary.tsv
```

The per-variant QC table and final rsid extract list are written alongside it:

```text
chr1_22_merged_hq.variant_qc.tsv
chr1_22_merged_hq.extract.txt
```

Finally, a second dsub job extracts the passing variants from the raw direct
bfile and writes:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq.{bed,bim,fam}
```

The final high-quality bfile keeps the same samples as the raw direct bfile;
only variants are filtered. The worker also computes sample missingness over
the final variant set for all samples and for the AoU EUR keep-list:

```text
chr1_22_merged_hq.sample_missingness_all.smiss
chr1_22_merged_hq.sample_missingness_eur.smiss
chr1_22_merged_hq.sample_missingness_summary.tsv
```

Observed Step 4 accounting from the default cdrv9/v9 run:

| Metric | Value |
|---|---:|
| AoU EUR keep-list samples | 302,714 |
| Raw direct-bfile variants | 501,452 |
| Raw direct-bfile samples | 535,662 |
| Final high-quality direct variants | 498,902 |
| Final high-quality direct samples | 535,662 |
| All-sample mean missingness | 0.0004140998315 |
| EUR-sample mean missingness | 0.0003748639908 |
| Maximum sample missingness | 0.0725814 |

Ordered v9 high-quality direct-SNP filter counts:

| Filter | Dropped SNPs | Remaining SNPs |
|---|---:|---:|
| Original UKBB direct-SNP rsid list | 0 | 501,477 |
| Present in the AoU direct bfile | 25 | 501,452 |
| Liftover ALT frequency available and alleles match | 0 | 501,452 |
| ALT-frequency difference <= 0.04 | 609 | 500,843 |
| AoU EUR MAF >= 0.007 | 6 | 500,837 |
| AoU EUR missingness <= 0.05 | 1,935 | 498,902 |

### Step 5 — ADMIXTURE K=6 projection

ADMIXTURE projection starts from the high-quality direct bfile:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq
```

`admixture_prep.sh` downloads two public inputs if they are not already cached
locally under gitignored paths: the ADMIXTURE Linux binary and the K=6 global
reference allele-frequency TSV from `public-statgen`. It stages those files to
the workspace bucket for Batch workers.

The prep worker then:

```text
1. computes all-sample variant missingness on the HQ direct bfile
2. keeps SNPs with all-sample missingness <= 0.05
3. intersects those SNPs with the K=6 reference SNPs
4. extracts the overlap into an ADMIXTURE-specific bfile
5. aligns PLINK BIM alleles to the reference frequencies
6. writes the aligned bfile plus ref_aligned.P
```

The all-sample missingness filter is specific to ADMIXTURE projection. It does
not modify the REGENIE Step-1 HQ direct bfile on disk.

Outputs from prep are written to:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/scrap/aou_admixture_aligned.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/scrap/ref_aligned.P
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/scrap/admixture_prep_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/scrap/admixture_align_log.txt
```

`admixture_split_batches.sh` splits the aligned bfile into batches of 20,000
people and copies `ref_aligned.P` to each batch as `batch_NNN.6.P.in`, which is
the filename ADMIXTURE expects in projection mode.

`admixture_run_projection.sh` submits one dsub task per missing batch Q file,
runs:

```text
admixture -j$(nproc) -P batch_NNN.bed 6
```

and concatenates the per-batch Q files into:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/aou_admixture_k6.tsv
```

The final TSV columns are:

```text
FID IID European East_Asian American African South_Asian Oceanian
```

Observed Step 5 accounting from the default cdrv9/v9 run:

| Metric | Value |
|---|---:|
| Source high-quality direct variants | 498,902 |
| Source samples | 535,662 |
| K=6 reference variants | 135,020 |
| Variants passing all-sample missingness <= 0.05 | 498,721 |
| Variants dropped by all-sample missingness | 181 |
| Variants after reference intersection | 134,748 |
| Variants dropped by allele alignment | 0 |
| Final aligned ADMIXTURE variants | 134,748 |
| Final aligned samples | 535,662 |
| `ref_aligned.P` rows | 134,748 |
| Projection batches | 27 |
| Full batch size | 20,000 samples |
| Final batch size | 15,662 samples |
| Projection Q rows | 535,662 |
| Final ADMIXTURE TSV samples | 535,662 |

### Step 6 — AoU-vs-ours ancestry comparison and European classifier

`compare_aou_ancestry.sh` always writes the European keep-list from this
pipeline's ADMIXTURE K=6 fractions. When an AoU RYE admixture-fraction file is
available, it also compares this pipeline's fractions to AoU's v8 auxiliary
ancestry outputs:

```text
aux/ancestry/echo_v4_r2.ancestry_preds.tsv
aux/admixture_estimates/aou_admixture_estimates_rye_v8.Q
```

AoU's RYE fractions include `mid` for Middle Eastern ancestry. Our K=6
projection uses the public-statgen global model and has `Oceanian`
instead of `mid`, so Step 6 treats MID as an AoU-specific component rather
than forcing a one-to-one mapping.

The European classifier uses fixed ancestry-fraction thresholds:

```text
European >= 0.8
African <= 0.1
American <= 0.1
East_Asian <= 0.1
Oceanian <= 0.1
South_Asian: no cap
```

The PLINK keep-list for downstream EUR-only analyses is written to:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/europeans/classified_european_iids.txt
```

The AoU-vs-ours comparison directory is:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/statgen/aou_vs_ours/
```

Full AoU RYE comparison summary outputs include:

```text
aou_vs_ours.summary.tsv
component_pair_metrics.tsv
component_correlation_matrix.tsv
fraction_distribution_summary.tsv
aou_pred_counts.tsv
aou_pred_vs_ours_european_summary.tsv
european_set_overlap_summary.tsv
european_set_group_ancestry_summary.tsv
discordant_set_component_summary.tsv
european_set_group_counts_by_aou_pred.tsv
aou_mid_threshold_summary.tsv
aou_mid_threshold_component_summary.tsv
aou_mid_bin_counts.tsv
aou_mid_bin_component_summary.tsv
```

When AoU RYE fractions are unavailable, Step 6 writes the European keep-list
plus classification-only summary tables:

```text
aou_vs_ours.summary.tsv
aou_pred_counts.tsv
aou_pred_vs_ours_european_summary.tsv
european_set_overlap_summary.tsv
```

The full comparison also writes plots under `aou_vs_ours/plots/`, including
component scatter plots, ancestry-fraction distributions, a hard-call vs
dominant-component heatmap, European set-overlap counts, discordant European
call composition plots, and AoU MID-threshold composition plots.

#### v9 default-run summary: AoU hard ancestry calls vs this pipeline's ADMIXTURE classifier

For v9, the AoU hard ancestry prediction file was available, but a v9 AoU RYE
admixture-fraction file was not found. Therefore this is a classification-only
comparison against AoU hard calls, not a full fraction-correlation validation.

The cdrv9 run included all 535,662 WGS samples present in the v9 ACAF files,
the AoU hard ancestry prediction file, and this pipeline's ADMIXTURE K=6
projection. The European classifier used the same fixed thresholds documented
above.

Key v9 European-set overlap results:

| Metric | Count |
|---|---:|
| AoU hard-predicted EUR | 302,714 |
| Classified European by this pipeline | 303,903 |
| European by both methods | 298,219 |
| AoU EUR only | 4,495 |
| This pipeline EUR only | 5,684 |
| Neither European | 227,264 |
| European union | 308,398 |

The v9 overlap was high: Jaccard index `0.966994`. Treating AoU's hard EUR
call as a reference label for validation only, this pipeline's European
classifier had precision `0.981297` and recall `0.985151`.

AoU hard-call counts in the v9 ancestry prediction file:

| AoU ancestry_pred | Count |
|---|---:|
| afr | 100,800 |
| amr | 107,928 |
| eas | 15,893 |
| eur | 302,714 |
| mid | 2,151 |
| sas | 6,176 |
| missing | 0 |

#### v8-only AoU RYE fraction-correlation validation

For v9, AoU provides an updated ancestry prediction file, but we did not find a
v9 equivalent of the AoU RYE admixture-fraction file. Therefore the full
AoU-vs-ours fraction comparison described here is a v8 validation analysis, not
a v9 pipeline requirement.

A historical v8 comparison used the AoU RYE fraction file and this pipeline's
ADMIXTURE K=6 projection. The European classifier used:

```text
European >= 0.8
African <= 0.1
American <= 0.1
East_Asian <= 0.1
Oceanian <= 0.1
South_Asian: no cap
```

Component-level correlations between AoU RYE fractions and this pipeline's
ADMIXTURE fractions were high for the directly comparable components:

| Component | Pearson r |
|---|---:|
| European | 0.9105 |
| East Asian | 0.9961 |
| American | 0.9495 |
| African | 0.9977 |
| South Asian | 0.9627 |

The European component has the lowest direct correlation because AoU RYE
includes a Middle Eastern (`mid`) component, while this pipeline's
public-statgen K=6 projection uses an Oceanian component instead. Step 6
therefore treats AoU MID as AoU-specific rather than forcing a one-to-one
component mapping.

### Step 7 — KING kinship and close relationship classification

Kinship starts from the high-quality direct bfile:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq
```

`subset_kinship_snps.sh` downloads the public UK Biobank SNP QC file
`ukb_snp_qc.txt` if it is not already cached locally, keeps SNPs where
`in_Relatedness == 1`, intersects those rsids with the HQ direct bfile, and
then submits a dsub worker to compute all-sample variant missingness on that
intersection. The final KING SNP list keeps variants with missingness strictly
less than `KINSHIP_MISSING_MAX` (default: `0.01`).

This is a **variant missingness** filter, not a sample missingness filter. No
`--mind`/sample filtering is applied before KING; the KING run uses all samples
present in `direct_bfile_hq/chr1_22_merged_hq`. The only sample-level output at
this stage is the Step 4 missingness report, which is diagnostic unless a later
pipeline step explicitly consumes it.

The subset stage writes:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/ukbb_relatedness_snps_in_hq_direct_geno_lt_threshold.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/kinship_snp_subset_all_sample_missingness.vmiss
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/kinship_snp_missingness_threshold_counts.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/kinship_snp_subset_summary.tsv
```

The summary reports the key set-operation counts:

```text
ukb_relatedness_snps
n_intersection_hq_direct
n_intersection_and_missing_lt_0.01
```

The threshold-count table reports how many SNPs would pass common missingness
cutoffs (`0.05`, `0.04`, `0.03`, `0.02`, `0.01`) from the same `.vmiss` file.
If the `.vmiss` file already exists and only `KINSHIP_MISSING_MAX` changes,
`subset_kinship_snps.sh` re-filters locally instead of re-submitting the dsub
job and re-localizing the large HQ direct bfile.

Observed Step 7a SNP-subset accounting from the default cdrv9/v9 run:

| Metric | Value |
|---|---:|
| Source high-quality direct variants | 498,902 |
| Source samples | 535,662 |
| UKBB `in_Relatedness` SNPs | 93,511 |
| Intersection with HQ direct bfile | 85,222 |
| Variants with all-sample missingness measured | 85,222 |
| Final KING SNPs after all-sample missingness <1% | 84,510 |

Common missingness-threshold counts from the same v9 `.vmiss` file:

| Threshold | Passing SNPs | Failing SNPs |
|---|---:|---:|
| missingness < 0.05 | 85,184 | 38 |
| missingness < 0.04 | 85,123 | 99 |
| missingness < 0.03 | 85,041 | 181 |
| missingness < 0.02 | 84,878 | 344 |
| missingness < 0.01 | 84,510 | 712 |

`get_genotypes.sh` proceeds directly from the SNP subset to the KING run. KING
is the largest single compute step in the pipeline, but it is idempotent: once
matching KING outputs exist, later runs skip it. To tighten the all-sample SNP
missingness threshold before KING, override `KINSHIP_MISSING_MAX` before
launching the pipeline:

```bash
# Override the default all-sample SNP missingness threshold.
KINSHIP_MISSING_MAX=0.02 bash get_genotypes.sh 2>&1
```

`run_king_kinship.sh` runs:

```text
plink2 --make-king-table --king-table-filter 0.035
```

on the HQ direct bfile with the final SNP extract list. It does not
materialize a separate kinship-only bfile; the exact SNP set is captured by
`ukbb_relatedness_snps_in_hq_direct_geno_lt_threshold.txt` and the parameter
files. KING writes:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/aou_hq_direct_rel.kin0
```

`kinship_qc.sh` compares our KING kinship coefficients against AoU's provided
KING-style relatedness table for the selected release:

```text
v9/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv
v8/wgs/short_read/snpindel/aux/relatedness/samples_relatedness.tsv
```

AoU's table provides the pairwise kinship coefficient, so this comparison is
kinship-only; IBS0 diagnostics are available only from our PLINK2/KING output.
Summary tables and plots are written under:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/qc/
```

`classify_relations.sh` then applies fixed KING kinship/IBS0 thresholds on our
KING table:

```text
identical:    kinship >= 0.3535 and IBS0 < 0.0012
parent_child: 0.1767 <= kinship < 0.3535 and IBS0 < 0.0012
sibling:      0.1767 <= kinship < 0.3535 and IBS0 >= 0.0012
```

The AoU version deliberately does not apply a birth-year/month sibling age-gap
filter because no portable AoU phenotype dependency is part of this genotype
pipeline.

Observed Step 7b-7d accounting from the default cdrv9/v9 run with
`KINSHIP_MISSING_MAX=0.01` and `KING_TABLE_FILTER=0.035`:

```text
UKBB in_Relatedness SNPs:                         93,511
Intersection with HQ direct bfile:                85,222
Final KING SNPs after all-sample missingness <1%: 84,510

KING pairs reported at kinship >=0.035:          116,490
AoU provided relatedness pairs:                   55,907
Overlapping pair set:                             55,744
Ours-only pair set:                               60,746
AoU-only pair set:                                   163
Pearson r vs AoU kinship:                         0.9896
Mean absolute kinship difference:                 0.00900
Median absolute kinship difference:               0.00599

Close relationships, kinship >=0.1767:            37,342
  sibling:                                        11,119
  parent_child:                                   22,248
  identical/twin/duplicate:                        3,975

Sibling relationship samples:
  total unique sibling individuals:                19,918
  pairs with both samples classified European:       4,651
  unique classified-European individuals in those pairs: 8,642
  classified-European individuals in any sibling pair:   8,671
  mixed sibling pairs with exactly one classified-European sample: 37
  classified-European individuals appearing only in mixed pairs:  29

Classified-European sibling-family components:
  families, defined as connected components of sibling pairs: 4,200
  exact family-size cells with 1-20 sibling individuals: suppressed
  non-suppressed exact family-size cells: 0
```

The two classified-European sibling-individual counts differ because one is
restricted to sibling pairs where both samples pass the European classifier,
whereas the other includes classified-European samples in mixed pairs where
only one sibling passes the classifier. In cdrv9, those 37 mixed sibling pairs
are threshold-edge cases: the non-European-classified sibling usually has an
`American` or `African` fraction just above the `0.1` cap, or a European
fraction just below `0.8`. The classifier is an operational GWAS keep-list
rule, not a family-level ancestry label.

These numbers are run accounting, not hardcoded expectations. New AoU releases
or changed thresholds should be summarized from the files in
`sbayesrc_genotypes/kinship/`.

### Step 8 — Select unrelated European IIDs for PCA fitting

`select_pca_europeans.sh` selects the European subset used to fit ancestry
PCA. It starts from this pipeline's European keep-list:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/europeans/classified_european_iids.txt
```

It then identifies European samples involved in `sibling` or `identical`
relationships from:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/close_relations.csv
```

Those samples are the seed exclusion set. The script expands that seed set to
include everyone directly related to those seeds at a third-degree kinship
threshold:

```text
KINSHIP >= 0.0441941
```

using:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/kinship/aou_hq_direct_rel.kin0
```

The candidate PCA set is:

```text
classified Europeans - expanded sibling/identical-relative exclusion set
```

Finally, the script runs PLINK2's maximal unrelated selector on the candidate
set:

```text
plink2 --psam candidate_pca_europeans.psam \
       --king-cutoff-table aou_hq_direct_rel.kin0 0.0441941
```

This final cutoff guarantees that no retained PCA-fitting pair has KING
kinship at or above the third-degree threshold. The step writes:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/fit_pca_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/select_pca_europeans.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_eur_log.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/seed_sibling_identical_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/expanded_sibling_identical_relatives_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/candidate_pca_europeans_iids.txt
```

The summary reports the number of classified Europeans, seed samples,
expanded exclusions, candidate PCA Europeans before the final KING cutoff,
samples removed by the final cutoff, and final `fit_pca_iids`.

Observed Step 8 accounting from the default cdrv9/v9 run with
`PCA_KINSHIP_THRESHOLD=0.0441941` and
`PCA_SEED_RELATIONSHIPS=sibling,identical`:

```text
Classified Europeans:                         303,903
Seed sibling/identical relationship rows:      15,094
Samples in seed relationship rows:             26,076
European seed samples:                         12,103
KING graph edges at kinship >=0.0441941:       87,885
KING graph nodes:                             116,528

Expanded exclusions total:                     13,971
Expanded exclusions in Europeans:              13,866
Expanded exclusions not in Europeans:             105
Candidate PCA Europeans before KING cutoff:   290,037
Removed by final KING cutoff:                  16,321
Final PCA fitting IIDs:                       273,716

Verification non-European retained:                 0
Verification expanded exclusions retained:          0
Verification related pairs retained:                0
```

These are accounting numbers from the cdrv9/v9 run, not hardcoded
expectations. New AoU releases or changed thresholds should be summarized from
`pca_eur/select_pca_europeans.summary.tsv`.

### Step 9 — QC SNPs for PCA

`pca_snp_qc.sh` builds the bfile used to fit European ancestry PCs. It starts
from the high-quality direct bfile:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq
```

and keeps only the final PCA-fitting IIDs from Step 8:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/fit_pca_iids.txt
```

The PCA-specific filters are applied in this order:

```text
1. keep final PCA-fitting IIDs
2. abs(AoU EUR ALT frequency - SBayesRC ALT frequency) <= 0.03
3. MAF >= 0.01 in the PCA-fitting samples
4. variant missingness <= 0.01
5. sample missingness <= 0.01
6. exclude long-range LD regions, hg38
7. LD prune with --indep-pairwise 1000 80 0.1
```

The tighter ALT-frequency agreement filter reuses Step 4's
`chr1_22_merged_hq.variant_qc.tsv`, which already contains
`abs_alt_freq_diff` after aligning AoU and SBayesRC hg38 alleles. This avoids
recomputing frequency concordance and keeps the Step 9 accounting directly
tied to the Step 4 QC table.

The long-range LD regions are hg38 regions downloaded from the public `plinkQC`
resource by the orchestrator and staged to the private Batch worker.

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_ready.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_snp_qc.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_snp_qc.filter_steps.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_afdiff_0.03.extract.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_ld_prune.prune.in
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_ld_prune.prune.out
```

`pca_snp_qc.filter_steps.tsv` reports, for every sequential filter, the input
SNP/sample count, dropped SNP/sample count, and remaining SNP/sample count.
`pca_snp_qc.summary.tsv` records the thresholds and final `pca_ready` counts.

Observed Step 9 accounting from the default cdrv9/v9 run:

| Filter | Dropped SNPs | Dropped samples | Remaining SNPs | Remaining samples |
|---|---:|---:|---:|---:|
| Source high-quality direct bfile | 0 | 0 | 498,902 | 535,662 |
| Keep Step 8 PCA-fitting IIDs | 0 | 261,946 | 498,902 | 273,716 |
| ALT-frequency difference <= 0.03 | 536 | 0 | 498,366 | 273,716 |
| MAF >= 0.01 | 3,482 | 0 | 494,884 | 273,716 |
| Variant missingness <= 0.01 | 3,561 | 0 | 491,323 | 273,716 |
| Sample missingness <= 0.01 | 0 | 135 | 491,323 | 273,581 |
| Exclude hg38 long-range LD regions | 8,029 | 0 | 483,294 | 273,581 |
| LD prune `1000 80 0.1` | 348,401 | 0 | 134,893 | 273,581 |

These numbers are run accounting, not hardcoded expectations. Recompute them
from `pca_snp_qc.filter_steps.tsv` after changing thresholds or moving to a
new AoU release.

### Step 10 — Fit PCA and project all samples

`fit_project_pca.sh` fits European ancestry principal components and projects
them onto the full high-quality direct bfile. The fit set is the Step 9
PCA-ready bfile:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_ready
```

That bfile contains the unrelated European PCA-fitting samples from Step 8
after the Step 9 SNP/sample QC and LD pruning. PCA is fit with PLINK2's
approximate PCA algorithm:

```text
plink2 --bfile pca_ready --pca allele-wts 20 approx --seed 0
```

The projection target is the all-sample high-quality direct bfile:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/direct_bfile_hq/chr1_22_merged_hq
```

The worker computes allele counts in the PCA fit set, verifies that every PCA
SNP is present in the projection bfile, parses the `A1` and PC-weight columns
from `aou_pcs.eigenvec.allele`, and projects all samples with PLINK2
`--score` using `no-mean-imputation` and `variance-standardize`.

PLINK2 writes `aou_pcs.eigenvec.allele` with one score row per allele. For
biallelic SNPs, that means two allele-weight rows per SNP; the pipeline
therefore validates the number of unique SNP IDs in that file, not the raw row
count.

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/aou_pcs.eigenval
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/aou_pcs.eigenvec
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/aou_pcs.eigenvec.allele
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/pca_eur_counts.acount
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/aou_projected.sscore
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/fit_project_pca.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/pca_eur/fit_project_pca.params.tsv
```

The summary reports fitted variants, PCA-fit samples, projected samples, PCs,
seed, score-file columns, and PCA SNPs missing from the projection bfile.
Missing projection SNPs should be zero because Step 9 starts from
`direct_bfile_hq`.

Observed Step 10 accounting from the default cdrv9/v9 run:

| Metric | Value |
|---|---:|
| PCA SNPs used for fitting/projection | 134,893 |
| PCA-fitting samples | 273,581 |
| Allele-weight rows | 269,786 |
| Unique allele-weight SNP IDs | 134,893 |
| Projected samples | 535,662 |
| Missing PCA SNPs in projection bfile | 0 |
| PCs | 20 |

Eigenvalues from the cdrv9/v9 run:

```text
539.662
110.526
87.4797
37.8953
30.6338
17.8471
13.1962
10.1802
9.91557
9.02512
8.94123
8.28208
8.17351
8.12065
8.08675
8.06028
8.02625
8.01715
7.99255
7.97727
```

### Step 11 — Sex covariate and sex/ploidy QC

`get_genetic_sex.sh` builds the binary sex covariate used by the height GWAS
example. It queries sex at birth from the AoU CDR `person` table, joins WGS
sex-ploidy metrics from AoU genomic QC, and keeps only samples with confident
binary sex:

```text
sex at birth is Female or Male
AND WGS sex ploidy is XX or XY
AND sex at birth and WGS ploidy are concordant
```

The output covariate uses FIDs from `direct_bfile_hq/chr1_22_merged_hq.fam`
and `sex_01` coding `0=female, 1=male`. The height setup step centers this as
`sex_c = sex_01 - 0.5`, so female is `-0.5` and male is `0.5`.

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/genetic_sex/sex_covar.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/genetic_sex/genetic_sex_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/genetic_sex/sex_ploidy_crosstab.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/genetic_sex/genetic_sex_log.txt
```

Observed Step 11 accounting from the default cdrv9/v9 run:

| Metric | Value |
|---|---:|
| Sample universe | 535,662 |
| Sex-at-birth query rows | 535,662 |
| Genomic metrics rows | 535,662 |
| Binary sex at birth | 530,385 |
| Missing sex-at-birth rows | 0 |
| Missing/non-binary sex at birth | 5,277 |
| Missing genomic metrics rows | 0 |
| Concordant confident binary samples | 527,995 |
| Confident female | 324,194 |
| Confident male | 203,801 |
| Binary sex/ploidy discordances | 493 |
| Non-canonical WGS sex ploidy | 1,896 |
| Confident sample percent | 98.568687% |

The cdrv9 run required concordance between sex at birth and WGS sex ploidy.

### Step 12 — Sample-QC exclusions for anomalous identical components

`build_identical_component_sample_qc.sh` builds a conservative GWAS exclusion
list from the close-relationship calls generated in Step 7. It treats
`relationship == identical` pairs as graph edges and computes connected
components. Ordinary monozygotic twin relationships are expected to appear as
simple size-2 components; components of size three or larger are not plausible
ordinary twin pairs and are excluded from downstream GWAS sample sets.

The default rule is:

```text
exclude all samples in identical-genotype connected components with size >= 3
```

The threshold can be changed with:

```bash
IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE=3
```

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/identical_components.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/identical_component_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/exclude_identical_component_size_ge3_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/identical_component_sample_qc.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/sample_qc/identical_component_sample_qc.log
```

Observed Step 12 accounting from the default cdrv9/v9 run:

| Metric | Value |
|---|---:|
| Close relationship rows | 37,342 |
| Identical KING pairs | 3,975 |
| Unique samples in identical graph | 6,434 |
| Identical components | 3,158 |
| Size-2 components | 3,088 |
| Size-3 components | 58 |
| Components size >=3 | 70 |
| Samples excluded | 258 |
| Complete-clique components | 3,158 |
| Non-complete components | 0 |
| Excluded samples missing from `.fam` | 0 |

Additional detailed component-size cells are suppressed in the README.

### Step 13 — Final GWAS genotype inputs

`make_gwas_genotype_inputs.sh` builds the final genotype files consumed by the
height REGENIE example. This step is separate from the earlier
`direct_bfile_hq` build because the final GWAS inputs use stricter,
purpose-specific QC:

```text
GWAS Step 1 source:
  direct_bfile_hq/chr1_22_merged_hq

GWAS Step 1 filters:
  variant missingness <= 0.01 in our classified European samples
  abs(fit_pca ALT frequency - SBayesRC/snp.info ALT frequency) <= 0.03
  MAF >= 0.007 in fit_pca_iids

GWAS Step 2 source:
  wgs_pfiles/chr{1..22}

GWAS Step 2 filters:
  variant missingness <= 0.03 in our classified European samples
  abs(fit_pca ALT frequency - SBayesRC/snp.info ALT frequency) <= 0.04
  MAF >= 0.007 in fit_pca_iids
```

The sample panels are intentionally different:

```text
classified European samples:
  europeans/classified_european_iids.txt
  used only for the genotype-missingness filter

fit_pca_iids:
  pca_eur/fit_pca_iids.txt
  used for ALT-frequency concordance and MAF
```

The SBayesRC/snp.info ALT frequency comes from
`data/support/sbayesrc_liftover_results.csv`. The filter helper matches
`A1_hg38`/`A2_hg38` to each AoU REF/ALT pair, converts `A1Freq` to the AoU ALT
allele, and drops a variant if the liftover row is missing or the alleles do
not match unambiguously.

The step first computes metrics with PLINK2:

```text
plink2 --bfile direct_bfile_hq/chr1_22_merged_hq \
  --keep pca_eur/fit_pca_iids.txt --freq counts
plink2 --bfile direct_bfile_hq/chr1_22_merged_hq \
  --keep europeans/classified_european_iids.txt --missing variant-only

plink2 --pfile wgs_pfiles/chrN \
  --keep pca_eur/fit_pca_iids.txt --freq counts
plink2 --pfile wgs_pfiles/chrN \
  --keep europeans/classified_european_iids.txt --missing variant-only
```

`build_gwas_genotype_filters.py` then applies the ordered filters, writes the
extract lists, writes per-variant QC tables, and writes a combined allele
frequency file for the passing 7M WGS variants. The MAF threshold replaces a
separate MAC threshold; the summary still reports the minimum fit-pca minor
allele count among passing variants as a sanity check.

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/gwas_genotype_qc.params.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/gwas_genotype_qc.summary.tsv

${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1.variant_qc.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1.fit_pca_alt_freqs.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/gwas_step1_direct.filter_steps.tsv

${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs_pfiles/chr{1..22}.{pgen,pvar,psam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs_pfiles/chr{1..22}.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/qc/chr{1..22}.variant_qc.tsv.gz
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/extracts/chr{1..22}.extract.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/fit_pca_af/chr{1..22}.fit_pca_alt_freqs_passing.tsv.gz
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/fit_pca_af/gwas_step2_fit_pca_alt_freqs_passing.tsv.gz
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs/gwas_step2_wgs.filter_steps.tsv
```

Observed Step 13 accounting from the completed cdrv9 run:

| Metric | Step 1 direct | Step 2 WGS |
|---|---:|---:|
| Source variants | 498,902 | 7,350,153 |
| Dropped: liftover missing / allele mismatch | 0 | 0 |
| Dropped: ALT-frequency difference above threshold | 500 | 23,135 |
| Dropped: fit-pca MAF < 0.007 | 0 | 232 |
| Dropped: classified-European missingness above threshold | 3,586 | 74,453 |
| Final variants | 494,816 | 7,252,333 |
| Final samples | 535,662 | 535,662 |
| Minimum fit-pca MAC among retained variants | 3,835 | 3,773 |
| Minimum fit-pca MAF among retained variants | 0.00700610 | 0.00700029 |
| Maximum classified-European missingness among retained variants | 0.0099999 | 0.0299997 |
| Maximum retained absolute ALT-frequency difference | 0.02999942 | 0.03999878 |

The drop counts are ordered/conditional: a variant removed by an earlier rule is
not counted again by later rules.

### Optional — HapMap3 HQ bfile

`make_hapmap3_bfile_hq.sh` is an optional post-pipeline helper. It builds a
European-sample HapMap3 bfile from the WGS pfiles using the existing Step 13
WGS allele-frequency and missingness metrics; it does not recompute those
metrics and it does not use the direct bfile as the variant source.

Run after Step 13 has completed:

```bash
bash make_hapmap3_bfile_hq.sh
```

Input HapMap3 rsids are tracked at the repo root:

```text
hapmap3_rsids.txt
```

Expected file properties:

```text
rows = 1,154,522
duplicates = 0
sha256 = 508e1b1739484b52af24b51f58aea833cf01f186b17696f0683ecd2abc687087
```

Variant filters are applied to HapMap3 rsids present in
`wgs_pfiles/chr{1..22}`:

```text
liftover allele/frequency available
alleles match liftover
abs(fit_pca ALT frequency - SBayesRC/snp.info ALT frequency) <= 0.03
MAF >= 0.007 in fit_pca_iids
variant missingness <= 0.01 in classified European samples
```

The final sample set is exactly:

```text
europeans/classified_european_iids.txt
```

No height phenotype, sex-covariate, or sample-QC exclusion filters are applied
to this optional bfile.

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/hapmap3_bfile_hq/hapmap3_bfile_hq.{bed,bim,fam}
${WORKSPACE_BUCKET}/sbayesrc_genotypes/hapmap3_bfile_hq/hapmap3_bfile_hq.filter_summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/hapmap3_bfile_hq/hapmap3_bfile_hq.variant_qc.tsv.gz
${WORKSPACE_BUCKET}/sbayesrc_genotypes/hapmap3_bfile_hq/hapmap3_bfile_hq.params.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/hapmap3_bfile_hq/hapmap3_bfile_hq.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/hapmap3_bfile_hq/hapmap3_bfile_hq.sample_missingness_eur.smiss
```

Observed optional HapMap3 HQ accounting from the completed cdrv9 run:

| Metric | Value |
|---|---:|
| Requested HapMap3 rsids | 1,154,522 |
| HapMap3 rsids present in WGS metrics | 1,154,381 |
| HapMap3 rsids absent from WGS metrics | 141 |
| Classified-European samples in final bfile | 303,903 |
| Fit-PCA samples used for frequency/MAF filters | 273,716 |
| Final HapMap3 HQ variants | 1,140,557 |
| Final HapMap3 HQ samples | 303,903 |
| Final `.bed` size | 81 GB |

### Step 14 — Height GWAS input setup

`setup_height_gwas.sh` builds the phenotype, covariate, and keep files for a
height GWAS in all samples classified as European by this pipeline. It uses
AoU program-collected height measurements:

```text
measurement_concept_id = 3036277
measurement_source_concept_id = 903133
measurement_type_concept_id = 44818701
unit_concept_id = 8582
minimum height = 140 cm
```

The BigQuery result is exported to the workspace bucket, processed inside the
AoU environment, and not downloaded off-platform. For each participant, the
helper takes the median valid height and the mean age at height measurement.
It then intersects with:

```text
1. our European ancestry keep-list
2. Step 12 sample-QC exclusion list
3. Step 11 confident sex covariates
4. the high-quality direct-bfile .fam sample IDs
5. Step 10 projected PC scores
```

Covariates are:

```text
age_c                  = age_at_height - mean(age_at_height)
sex_c                  = sex_01 - 0.5
age_c_sex_c_inter      = age_c * sex_c
PC1_AVG ... PC10_AVG   = first 10 projected European PCs
```

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/training_iids.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/phen.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/covar.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/base_covar.txt
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/height_gwas.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_input/height_example/height_gwas_log.txt
```

Observed Step 14 accounting from the completed cdrv9 run:

| Metric | Value |
|---|---:|
| Our classified Europeans | 303,903 |
| Samples in sample-QC exclusion list | 258 |
| Europeans removed by sample QC | 96 |
| Europeans after sample QC | 303,807 |
| Confident sex covariate rows | 527,995 |
| Height query rows after source/min-height filters | 507,036 |
| Projected PC rows | 535,662 |
| Europeans missing a height row | 54,336 |
| Height candidates missing confident sex | 3,578 |
| Height+sex candidates missing fam row | 0 |
| Height+sex+fam candidates missing PCs | 0 |
| Final GWAS samples | 245,893 |
| GWAS female | 146,124 |
| GWAS male | 99,769 |
| Mean height | 169.193404 cm |
| Median height | 168.4 cm |
| Mean age at height measurement | 55.7595201 years |
| PCs included | 10 |

### Step 15 — Height GWAS with REGENIE

`run_continuous_regenie_gwas.sh` runs the optional continuous-trait REGENIE
height example. It is gated by `RUN_HEIGHT_GWAS=1` because it launches one
Step 1 Batch job plus a 22-task Step 2 Batch array. The default run applies
rank-inverse normal transformation (`--apply-rint`) and uses the Step 14
covariates.

Step 1 uses the final GWAS direct-SNP bfile from Step 13:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step1_direct/chr1_22_merged_gwas_step1
```

Step 2 tests the final Step 13 WGS pfiles:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/gwas_genotypes/step2_wgs_pfiles/chr{1..22}
```

The orchestrator stages a small REGENIE runtime bundle to the workspace bucket
from the Workbench-provided REGENIE binary. The bundle includes the MKL shared
libraries needed on the default Ubuntu Batch worker. For AoU WGS pfiles, the
Step 2 worker also writes a local REGENIE-compatible `.psam` with `FID=0`
because AoU pfiles use a `#IID`-only psam header. It rewrites the localized
Step 1 prediction list so Step 2 reads the localized LOCO prediction file.
Neither change modifies the stored pfiles.

Outputs:

```text
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/step1/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/step2/chr{1..22}/
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/lightweight/chr{1..22}.height_example.regenie_lite.tsv.gz
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/lightweight/regenie_lite.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/regenie_gwas.summary.tsv
${WORKSPACE_BUCKET}/sbayesrc_genotypes/regenie_output/height_example/regenie_gwas.params.tsv
```

After Step 2 verification, the shared runner writes compact per-chromosome
association summaries with columns `rsid`, `allele1`, `a1freq`, `n`, `beta`,
`se`, and `log10p`. Set `REGENIE_MAKE_LIGHTWEIGHT_OUTPUTS=0` or pass
`--no-lightweight` to skip these compact files.

Observed Step 15 accounting from the completed cdrv9 run with
`RUN_HEIGHT_GWAS=1`:

| Metric | Value |
|---|---:|
| Step 1 analysis samples | 245,893 |
| Step 1 direct-SNP variants | 494,816 |
| Step 2 WGS chromosomes | 22 |
| Step 2 tested WGS variants | 7,252,333 |
| Lightweight association rows | 7,252,333 |
| Rank-inverse normal transform applied | yes |
| PCs included as covariates | 10 |

## Downstream workflows

These optional workflows build on the genotype preparation, sample selection
and ancestry analysis above. They have separate entry points and output
directories; `get_genotypes.sh` does not launch them automatically.

### `pan_aou_gwas/`: phenotype-wide association testing

[`pan_aou_gwas/`](pan_aou_gwas/README.md) builds many phenotypes from AoU survey
responses, physical measurements and ZIP3 socioeconomic context, then runs
PLINK 2 association tests using the shared HapMap3 HQ genotype set. It uses the
third-degree-unrelated European `fit_pca_iids` set, with the workflow's sex
covariate and identical-component exclusions. Phenotype manifests record
coding, sample restrictions and covariates; quantitative and binary traits are
pre-residualized before linear association testing.

This is useful for producing consistently constructed GWAS summary statistics
across many traits without rebuilding genotypes for each phenotype. Complete
the main pipeline's sample/PC inputs and the [optional HapMap3 HQ
bfile](#optional--hapmap3-hq-bfile) first. See its README for setup-only, smoke
test and full-run commands, phenotype definitions and analysis assumptions.

### `phasing_and_gnomix/`: haplotype ancestry and scoring inputs

[`phasing_and_gnomix/`](phasing_and_gnomix/README.md) selects all WGS samples
outside `classified_european_iids.txt`, together with a reproducible random
sample of 50,000 unrelated Europeans from `fit_pca_iids.txt`. It reads the
original AoU ACAF callset and phases the union of lifted pretrained Gnomix
features and `sbayesrc_hg38.csv` SNPs, retaining matched variants with at most
10% missing genotypes. Beagle performs population phasing; the pretrained
Gnomix model and Gnofix then run for every selected sample.

Unlike genome-wide ADMIXTURE fractions, local ancestry identifies the ancestry
assigned to each haplotype along a chromosome. This enables downstream analyses
such as applying an ancestry-specific predictor to segments of that ancestry.
The workflow prepares the required inputs; it does not fit predictors or
calculate polygenic scores. Outputs include Gnofix ancestry tracts, ancestry
fractions, window-level called-ancestry posteriors, and retained SBayesRC SNPs
with Gnofix phase corrections applied. All eight ancestry labels remain
separate, including EUR and WAS.

The workflow uses chromosome/shard directories in the workspace GCS bucket,
tracked seeds, pinned public references and a resumable Batch coordinator.
It validates a full chr22 pilot before submitting the remaining chromosomes.
The README documents the output layout, allele/window mappings, resource
settings and validation status. Participant data, run logs and private notes
remain outside the public source directory; no karyograms are generated.

### Other phenotype workflows

Additional downstream workflows for educational attainment, household income,
ETM cognitive task scores, and the final GradCPT/Flanker-enriched proxy
phenotype are documented in [ea_proxy.md](ea_proxy.md). They are separate from
the core genotype-preparation pipeline and assume `get_genotypes.sh` has
completed through the final REGENIE genotype-input step.

## Prerequisites

You are inside an AoU Verily Jupyter session terminal (the standard
interactive analysis environment for the All of Us Researcher Workbench).

Before cloning the repo in a brand-new workspace, add the current Controlled
Tier data collection from the Workbench catalog:

1. Open **Data from Catalog** with the plus button in the workspace.
2. Select **All of Us Controlled Tier**.
3. On **Select resource**, use the recommended active version (`cdrv9`) and
   choose **Select all resources**.
4. On **Review policies**, answer the Research Use Statement questions and
   acknowledge the permanent Controlled Tier workspace policies.
5. On **Review selection**, add all six selected resources at the workspace
   root. The selected resources should include the `v9-genomics-folder` object
   and the `vwb-aou-datasets-controlled-v9` bucket reference.

After the resources are added and the Jupyter environment is running, the
clean-workspace terminal workflow is:

```bash
git clone https://github.com/human-genomics/aou-sbayesrc-gwas.git
cd aou-sbayesrc-gwas
bash get_genotypes.sh
```

or, to include the optional height GWAS:

```bash
RUN_HEIGHT_GWAS=1 nohup bash get_genotypes.sh > logs/run_height_gwas.log 2>&1 &
```

`get_genotypes.sh` handles first-run setup that can safely be automated from
inside the AoU session: local output directories, Python dependencies, public
reference downloads, the Workbench workspace-bucket resource/mount when
possible, and the temporary BigQuery dataset used by CDR-querying steps. If
AoU or Workbench permissions prevent creating/mounting those managed
resources, the script fails before doing pipeline work and prints the missing
prerequisite.

The script selects the AoU release explicitly with `AOU_DATA_VERSION` (`v9` by
default). For v9 it sets `WORKSPACE_CDR` to `C2025Q4R6` internally, even if the
Workbench session exported an older default CDR. Set `AOU_STRICT_WORKSPACE_CDR=1`
to fail instead of overriding a mismatched `WORKSPACE_CDR`.

The session must provide:

- **Env vars** set automatically by the Workbench: `GOOGLE_PROJECT`,
  `CDR_STORAGE_PATH`. (`$WORKSPACE_BUCKET` is unreliable on some pods —
  `get_genotypes.sh` derives the writable bucket URI from the FUSE mount
  table instead, which is portable across users.)
- **Controlled-tier dataset bucket** FUSE-mounted (read-only) at
  `/home/jupyter/workspace/vwb-aou-datasets-controlled-v9/`.
- **Workspace bucket** FUSE-mounted (read-write) at
  `/home/jupyter/workspace/workspace-bucket/`. If it is absent,
  `get_genotypes.sh` attempts to create the Workbench-controlled
  `workspace-bucket` GCS resource and run `wb resource mount`, then verifies
  that the path is a real writable `gcsfuse` mount. It intentionally does not
  create this path as a local directory, because that would write large
  pipeline outputs to ephemeral disk instead of durable workspace storage.
- **`plink2`** preinstalled at `/opt/workbench-tools/binaries/bin/plink2`.
  The orchestrator stages this binary to the workspace bucket once, and each
  Batch worker `--input`s it back.
- **`regenie`** preinstalled at `/opt/workbench-tools/binaries/bin/regenie`
  for the optional height GWAS example.
- **`dsub`** preinstalled at `/opt/conda/envs/jupyter/bin/dsub`.
- **Internet access from the Jupyter pod** for public reference downloads.
  Batch workers use private networking, so the orchestrator downloads public
  ADMIXTURE inputs and stages them to the workspace bucket.
- **gcloud Application Default Credentials** (run `gcloud auth
  application-default login` once if `gcloud storage` or dsub return auth
  errors).
- **Pet service account** of the pod, used as the Batch worker identity
  (auto-detected via `gcloud config get-value account`).

Python dependencies for the local helper scripts are listed in
`requirements.txt` and auto-installed by `get_genotypes.sh`. `dsub` is provided
by the Workbench image. CDR-querying steps use `SBAYESRC_BQ_TMP_DATASET` or
their step-specific override if set; otherwise they reuse an existing dataset
or attempt to create `sbayesrc_tmp`.

## Usage

Smoke-test one chromosome of the WGS extraction step (one dsub task, ~15–25 min;
the direct-bfile step is skipped because it requires all 22 chromosomes):
```bash
SBAYESRC_TEST_CHROM=22 bash get_genotypes.sh 2>&1
```

Full all-22 run (recommended in the background):
```bash
nohup bash get_genotypes.sh > logs/run.log 2>&1 &
```

Run the optional height GWAS example as well:
```bash
RUN_HEIGHT_GWAS=1 nohup bash get_genotypes.sh > logs/run_height_gwas.log 2>&1 &
```

Build the optional HapMap3 HQ bfile after Step 13 has completed:
```bash
bash make_hapmap3_bfile_hq.sh
```

The default `bash get_genotypes.sh` command is designed to run the full
pipeline end to end. Expensive steps, including KING and the optional REGENIE
height GWAS, have parameter/count-based idempotency checks so later runs skip
matching outputs instead of re-submitting work.

The run logs to `logs/run_YYYYMMDD_HHMMSS.log` (timestamped) and tees through
to the foreground if attached. Each Batch worker's stdout/stderr is uploaded
to `${WORKSPACE_BUCKET}/sbayesrc_genotypes/logs/dsub/`.

Steps that query the AoU CDR use a temporary BigQuery dataset in the user's
workspace project. The pipeline uses `SBAYESRC_BQ_TMP_DATASET` if set;
otherwise it chooses an existing dataset in this order: `sbayesrc_tmp`,
`high_quality_cohort`, then the first dataset returned by `bq ls`. If no
dataset exists, CDR-querying setup steps attempt to create `sbayesrc_tmp`
through Workbench and then through `bq`. The AoU pet service account may not be
allowed to create new BigQuery datasets, so set `SBAYESRC_BQ_TMP_DATASET` to
an existing writable dataset if automatic creation is blocked or
auto-selection is not appropriate for your workspace.

## Portability

No user-specific identifiers are hardcoded. The pipeline derives:

- the **workspace bucket URI** from the FUSE mount table (not `$WORKSPACE_BUCKET`, which on some pods points to a non-existent cloned-bucket name);
- the **GCP project** from `$GOOGLE_PROJECT` (set by the Workbench);
- the **pet service account** from `gcloud config get-value account` (specific to each pod);
- the **VPC network/subnetwork** as `projects/$GOOGLE_PROJECT/global/networks/network` and `projects/$GOOGLE_PROJECT/regions/us-central1/subnetworks/subnetwork` — these are the canonical AoU project VPC names and are the same for every user;
- the **AoU controlled-tier path** as a stable gs:// URI shared across users.

Any AoU researcher with controlled-tier access should be able to clone the
repo and run it as-is.

## Idempotency

- Step 1 skips chromosomes whose local `chr{N}.{extract,idmap}.txt` already
  exist non-empty.
- Step 2 skips chromosomes whose `chr{N}.pgen` already exists at the
  destination (cheap gcsfuse metadata-only lookup via `test -f` on the
  workspace bucket mount — no API call needed for the skip check).
- Step 3 skips direct pfiles and the merged direct bfile when they already
  exist with the expected variant counts.
- Step 4 skips the high-quality direct bfile when its parameter file, summary,
  bfile, and sample-missingness summary already exist with matching thresholds
  and expected variant counts. The EUR frequency/missingness metrics are also
  skipped when their output line counts match the raw direct bfile.
- Step 5 skips ADMIXTURE prep, split, projection batches, and final concat
  independently when their parameter files and expected row counts match.
  Projection only submits missing or stale batch `.Q` files.
- Step 6 skips ancestry classification/comparison when its parameter file,
  required summary tables, required plots (for full AoU RYE comparison mode),
  and European keep-list already exist with matching input sizes and
  thresholds.
- Step 7 skips the kinship SNP subset, KING run, AoU comparison, and close
  relationship classifier independently when their parameter files, thresholds,
  and expected counts match.
- Step 8 skips the PCA European selector when its parameter file, summary, and
  `fit_pca_iids.txt` exist with matching input sizes, threshold, seed
  relationships, and output count.
- Step 9 skips PCA SNP QC when its parameter file, summary, filter-step table,
  and `pca_ready.{bed,bim,fam}` exist with matching input sizes, thresholds,
  high-LD-region checksum, and final counts.
- Step 10 skips PCA fitting/projection when its parameter file, summary,
  eigenvalue/eigenvector/allele-weight files, allele counts, and projected
  score file exist with matching input sizes, PC settings, and output counts.
- Step 11 skips the sex covariate/QC build when its parameter file, summary,
  crosstab, log, and `sex_covar.txt` exist with matching input sizes and
  concordance settings.
- Step 12 skips identical-component sample QC when its parameter file, summary,
  component tables, exclusion list, and log exist with matching input sizes
  and exclusion threshold.
- Step 13 skips final GWAS genotype construction when its parameter file,
  summary, extract lists, Step 1 bfile, and Step 2 per-chromosome pfiles exist
  with matching input sizes, thresholds, and expected variant counts. It also
  skips metric-scan dsub jobs independently when their output line counts match
  the source genotype files and sample keep-list sizes.
- Step 14 skips height GWAS input setup when its parameter file, summary,
  phenotype, covariate, keep, and log files exist with matching concept IDs,
  minimum-height threshold, PC count, input sizes, and sample-QC exclusion
  list.
- Step 15 is skipped unless `RUN_HEIGHT_GWAS=1`. When enabled, it skips
  REGENIE Step 1 if the prediction list, params, and summary match the desired
  run parameters. It then submits only Step 2 chromosomes missing matching
  per-chromosome summary/params files, and writes a final GWAS summary only
  after all selected chromosomes validate.

A re-run of `get_genotypes.sh` after a successful run should skip every step
and submit zero dsub tasks.

## Files

| File | Role |
|---|---|
| `get_genotypes.sh` | Orchestrator. Sets paths/env (including the dsub knobs) and calls each step script. |
| `store_sbayesrc_ids.py` | Step 1 — local SBayesRC ID + idmap generation. |
| `wgs_extract_variants.sh` | Step 2 — submits the dsub `--tasks` TSV (one row per missing chrom), waits, polls dstat until every task is terminal, verifies bucket outputs, writes combined summary. |
| `dsub_extract_worker.sh` | Per-Batch-worker script — the four-pass extract for one chromosome. |
| `normalize_pvar_alleles.py` | Canonical, testable reference for the pass-1.5 right-trim logic that's inlined as awk in the worker. Not used at runtime — kept for unit-testability and local debugging. |
| `prepare_direct_snps.py` | Step 3a — prepares per-chromosome direct-SNP extraction lists and missing-SNP reports. |
| `extract_direct_snps.sh` | Step 3b — extracts present direct SNPs from WGS pfiles. |
| `make_direct_bfile.sh` | Step 3c — submits/verifies the direct-bfile dsub merge. |
| `dsub_direct_bfile_worker.sh` | Per-Batch-worker script — merges direct pfiles into a sorted intermediate pgen, converts to `chr1_22_merged.{bed,bim,fam}`, and writes a summary. |
| `make_hq_direct_bfile.sh` | Step 4 — creates the AoU EUR keep-list, submits EUR frequency/missingness metrics, applies high-quality direct-SNP filters, builds the final filtered bfile, and verifies outputs. |
| `dsub_hq_direct_metrics_worker.sh` | Step 4 metrics worker — computes EUR allele frequencies and variant missingness on the raw direct bfile. |
| `filter_hq_direct_snps.py` | Step 4 local filter builder — joins EUR metrics to SBayesRC liftover frequencies and writes the ordered filter summary, variant QC table, and final extract list. |
| `dsub_hq_direct_bfile_worker.sh` | Step 4 bfile worker — extracts the passing variants, writes `chr1_22_merged_hq.{bed,bim,fam}`, and computes sample missingness over the final variant set. |
| `admixture_prep.sh` | Step 5a — stages ADMIXTURE inputs and submits the prep/alignment worker. |
| `dsub_admixture_prep_worker.sh` | Step 5a worker — applies all-sample `geno <= 0.05`, intersects reference SNPs, aligns alleles, and writes `aou_admixture_aligned` plus `ref_aligned.P`. |
| `admixture_align_alleles.py` | Step 5a reference helper — readable/testable Python version of the allele-alignment logic; the Batch worker uses an inline awk implementation to avoid depending on python inside the worker image. |
| `admixture_split_batches.sh` | Step 5b — submits the batch-splitting worker. |
| `dsub_admixture_split_worker.sh` | Step 5b worker — creates 20,000-person batch bfiles and per-batch `.P.in` files. |
| `admixture_run_projection.sh` | Step 5c — submits projection tasks and the final concat job. |
| `dsub_admixture_project_worker.sh` | Step 5c worker — runs ADMIXTURE `-P` for one batch and verifies Q row counts. |
| `dsub_admixture_concat_worker.sh` | Step 5c worker — concatenates batch Q files and writes `aou_admixture_k6.tsv`. |
| `compare_aou_ancestry.sh` | Step 6 — idempotent wrapper for the AoU-vs-ours ancestry comparison and European keep-list. |
| `compare_aou_ancestry.py` | Step 6 helper — joins AoU hard ancestry calls, AoU RYE fractions, and our ADMIXTURE fractions; writes aggregate tables and plots. |
| `classify_admixture_europeans.py` | Step 6 fallback helper — writes the European keep-list and hard-call overlap summaries when AoU RYE fractions are unavailable. |
| `subset_kinship_snps.sh` | Step 7a — intersects UKBB `in_Relatedness == 1` SNPs with the HQ direct bfile and applies all-sample missingness `< KINSHIP_MISSING_MAX`. |
| `dsub_kinship_subset_worker.sh` | Step 7a worker — computes variant missingness and writes the final KING SNP extract list plus count summary. |
| `run_king_kinship.sh` | Step 7b — submits/verifies the PLINK2 KING run from the HQ direct SNP subset. |
| `dsub_king_kinship_worker.sh` | Step 7b worker — runs `plink2 --make-king-table --king-table-filter 0.035`. |
| `kinship_qc.sh` | Step 7c — idempotent wrapper comparing our KING output to AoU's provided relatedness table. |
| `kinship_qc.py` | Step 7c helper — writes kinship comparison summaries, pair-level overlap data, scatter plots, and Bland-Altman plots. |
| `classify_relations.sh` | Step 7d — idempotent wrapper for close-relationship classification. |
| `classify_relations.py` | Step 7d helper — classifies sibling, parent_child, and identical pairs from KING kinship/IBS0 thresholds. |
| `select_pca_europeans.sh` | Step 8 — idempotent wrapper selecting unrelated European IIDs for fitting PCA. |
| `select_pca_europeans.py` | Step 8 helper — expands European sibling/identical seeds to their third-degree relatives, removes them, and runs PLINK2 `--king-cutoff-table` for the final unrelated set. |
| `pca_snp_qc.sh` | Step 9 — stages PCA SNP QC inputs and submits/verifies the Batch worker that builds `pca_ready.{bed,bim,fam}`. |
| `dsub_pca_snp_qc_worker.sh` | Step 9 worker — applies the PCA-fitting sample keep-list, tighter frequency concordance, MAF/geno/mind filters, long-range-LD exclusion, and LD pruning. |
| `fit_project_pca.sh` | Step 10 — stages PCA inputs and submits/verifies the Batch worker that fits European PCs and projects them to all samples. |
| `dsub_fit_project_pca_worker.sh` | Step 10 worker — fits PLINK2 PCA with allele weights, computes fit-set allele counts, verifies projection SNP coverage, and writes all-sample PC scores. |
| `get_genetic_sex.sh` | Step 11 — queries AoU sex at birth, joins WGS sex ploidy, and writes the confident binary sex covariate plus QC summaries. |
| `get_genetic_sex.py` | Step 11 helper — builds `sex_covar.txt`, the sex/ploidy crosstab, summary, and verification log. |
| `build_identical_component_sample_qc.sh` | Step 12 — builds the sample-QC exclusion list for identical-genotype components of size `IDENTICAL_COMPONENT_EXCLUDE_MIN_SIZE` or larger. |
| `build_identical_component_sample_qc.py` | Step 12 helper — computes identical-pair connected components, writes component tables, and writes the exclusion FID/IID file. |
| `make_gwas_genotype_inputs.sh` | Step 13 — computes GWAS-specific allele frequency/missingness metrics, builds final Step 1/Step 2 extract lists, and submits the final genotype-build workers. |
| `build_gwas_genotype_filters.py` | Step 13 helper — joins fit-pca allele counts, classified-EUR missingness, and SBayesRC/snp.info frequencies; writes filter summaries, variant QC tables, extract lists, and passing Step 2 fit-pca ALT frequencies. |
| `dsub_gwas_direct_metrics_worker.sh` | Step 13 worker — computes fit-pca allele counts and classified-EUR missingness for `direct_bfile_hq`. |
| `dsub_gwas_wgs_metrics_worker.sh` | Step 13 worker — computes fit-pca allele counts and classified-EUR missingness for one WGS pfile chromosome. |
| `dsub_gwas_step1_direct_worker.sh` | Step 13 worker — extracts the final REGENIE Step 1 bfile from `direct_bfile_hq`. |
| `dsub_gwas_step2_wgs_worker.sh` | Step 13 worker — extracts one final REGENIE Step 2 WGS pfile chromosome. |
| `make_hapmap3_bfile_hq.sh` | Optional helper — builds a classified-European HapMap3 bfile from WGS pfiles using existing Step 13 WGS QC metrics. |
| `filter_hapmap3_wgs_hq_snps.py` | Optional helper — filters tracked HapMap3 rsids by liftover allele match, fit-pca frequency/MAF, and classified-European missingness. |
| `dsub_hapmap3_wgs_extract_worker.sh` | Optional worker — extracts one HapMap3 HQ chromosome pfile from WGS pfiles and keeps classified-European samples. |
| `dsub_hapmap3_bfile_merge_worker.sh` | Optional worker — merges HapMap3 HQ chromosome pfiles and converts them to a single PLINK bfile. |
| `setup_height_gwas.sh` | Step 14 — queries program-collected AoU height, exports the result inside the workspace bucket, and builds REGENIE phenotype/covariate/keep files after sample-QC exclusions. |
| `setup_height_gwas.py` | Step 14 helper — intersects Europeans, sample-QC exclusions, height, confident sex, genotype IDs, and projected PCs; centers covariates; writes summaries and verification checks. |
| `run_continuous_regenie_gwas.sh` | Step 15 — optional continuous-trait REGENIE runner using the final Step 13 genotype inputs. |
| `make_lightweight_regenie_outputs.py` | Step 15 helper — converts full per-chromosome REGENIE outputs into compact `rsid/allele/frequency/effect/SE/log10p` TSVs. |
| `dsub_regenie_step1_worker.sh` | Step 15 worker — runs REGENIE Step 1, writes LOCO predictions, and verifies sample/variant counts. |
| `dsub_regenie_step2_worker.sh` | Step 15 worker — runs one REGENIE Step 2 chromosome, adapting AoU `#IID` psam headers and localized Step 1 prediction paths for REGENIE. |
| `requirements.txt` | Python dependencies for the local helper scripts. |
| `CLAUDE.md` | Project conventions, AoU platform notes, portability rules, dsub-from-Jupyter recipe. Gitignored — local developer reference. |
