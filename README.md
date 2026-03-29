# nxf-susie-coloc
Performs Finemapping and colocalization of GWAS results using the the provided LD reference dataset

## TODO

- [ ] Make locus_kb a Nextflow argument

## Running the Pipeline

```bash
nextflow run main.nf -profile eddiedev -resume -c run.config
```

## Development

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