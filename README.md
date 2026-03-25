# nxf-susie-coloc
Performs Finemapping and colocalization of GWAS results using the the provided LD reference dataset

## Development

**Singularity**

Example command with mount point:

```bash
singularity shell \
    -B /gpfs/igmmfs01/eddie/ISARIC4C/GTEx_Analysis_v10_QTLs/:/mnt/gtex \
    -B /home/olabayle/isaric/olivier/Covid19/data:/mnt/gwas_data \
    -B $PWD:/FinemapColoc \
    docker://olivierlabayle/nxf-susie-coloc:main
```