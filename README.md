# nxf-susie-coloc

Performs Finemapping and colocalization of GWAS results using the the provided LD reference dataset

## Pipeline Arguments

### Compurlsory Arguments

- `GWAS_OUTCOME_VAR (default: 0.16)`: Variance of the GWAS outcome.
- `GWAS_OUTCOME_TYPE (default: "cc")`: GWAS outcome type ("cc" is binary).
- `GTEX_FILES (default: "/gpfs/igmmfs01/eddie/ISARIC4C/GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations/*.parquet")`: Location of GTEX parquet files.
- `GTEX_TISSUES (default: ["Whole_Blood", "Lung", "Spleen", "Heart_Left_Ventricle", "Heart_Atrial_Appendage", "Artery_Coronary", "Artery_Aorta", "Minor_Salivary_Gland", "Esophagus_Mucosa", "Colon_Sigmoid", "Colon_Transverse", "Liver", "Adipose_Visceral_Omentum", "Adipose_Subcutaneous", "Kidney_Cortex"])`: GTEX tissues to colocalize.
- `GWAS_RESULTS defaukt: "/gpfs/igmmfs01/eddie/ISARIC4C/olivier/GenOMICC_Mortality/data/META_ANALYSIS.ALIVE_AT_ASSESSMENT.gwas.tsv"`: Location of GWAS summary statistics.
- `REFERENCE_PREFIX (default: "/gpfs/igmmfs01/eddie/ISARIC4C/olivier/data/kgp-merged-unrelated-or3/kgp.merged.unrelated")`: Location of reference panel in plink bed format.

### GWAS Clump Identification Options

- `GWAS_FP_DIRS (default: "results_backup/gwas_fp_results/*")`: If clumps have already been identified separately, the pipeline can start from there.
- `GWAS_MIN_SIG_CLUMP_SIZE (default: 7)`: Minimum number of variants in a clump to be considered significant.
- `GWAS_LEAD_PVALUE (default: 5e-8)`: Lead variant p-value.
- `GWAS_P2_PVALUE (default 5e-5)`: Secondary p-value used by plink2.
- `GWAS_R2_THRESHOLD (default: 0.2)`: Minimum R2 in the clump.
- `GWAS_CLUMP_KB (default: 500)`: Size of the clump.

### GTEX Options

- `GTEX_SAMPLE_SIZE = 940`: Sample size of GTEX.

### Finemapping Options

- `SUSIE_COVERAGE (default: 0.8)`: Susie coverage argument.
- `SUSIE_MAXIT (default: 1000)`: Susie maximum number of iterations.
- `LOCUS_KB (default: 250)`: finemapping locus size.

### Development Options

- `USE_SYSIMAGE (default: true)`: Whether to use the Julia sysimage, only set to false when developping.

## Running the Pipeline

Update the Nextflow profile and `run.config` file as desired, then run.

```bash
nextflow run main.nf -profile eddiedev -resume -c run.config
```

## Development Tips

**Singularity**

Example command with mount point:

```bash
SINGULARITY_DISABLE_CACHE=1 singularity shell \
    -B /gpfs/igmmfs01/eddie/ISARIC4C/GTEx_Analysis_v10_QTLs/:/mnt/gtex \
    -B /home/olabayle/isaric/olivier/data:/mnt/gwas_data \
    -B $PWD/src:/opt/FinemapColoc/src \
    -B $PWD/results:/mnt/results \
    --no-home \
    docker://olivierlabayle/nxf-susie-coloc:main
```

Example of command run:

```bash
JULIA_DEPOT_PATH=/tmp:\$JULIA_DEPOT_PATH julia --project=/opt/FinemapColoc \
    /opt/FinemapColoc/bin/run.jl prepare-gwas-results \
    /mnt/gwas_data/Covid19/covid_19_results_2026/meta_analysis_workdir/META_ANALYSIS.all.tsv /mnt/gwas_data/kgp-merged-unrelated-or3/kgp.merged.unrelated
```