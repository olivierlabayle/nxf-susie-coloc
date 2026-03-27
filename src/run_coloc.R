library(coloc)
library(susieR)
library(data.table)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

gwas_susie_results <- readRDS(args[1])
qtl_susie_results <- readRDS(args[2])
output_prefix <- args[3]

# gwas_susie_results <- readRDS("/gpfs/igmmfs01/eddie/ISARIC4C/olivier/nxf-susie-coloc/work/ac/fa11b2d1cd7a00ed4e001af39cc5d1/gwas_fp_results_chr17_67835900/GWAS.chr17.67835900.susie_results.rds")
# qtl_susie_results <- readRDS("GTEX_17_67835900_Lung_ENSG00000237854.4/GTEX_17_67835900_Lung_ENSG00000237854.4.susie_results.rds")
# output_prefix <- "GTEX_17_67835900_Lung_ENSG00000237854.4/GTEX_17_67835900_Lung_ENSG00000237854.4.coloc"

coloc_results = coloc.susie(gwas_susie_results, qtl_susie_results)

if (any(is.na(coloc_results))){
    print("Could not colocalize")
} else {
    print("Could colocalize.")
    saveRDS(coloc_results, file = paste0(output_prefix, ".rds"))
    fwrite(coloc_results$summary, paste0(output_prefix, ".tsv"), sep = "\t")
}
