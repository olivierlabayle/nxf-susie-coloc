library(data.table)
library(coloc)
library(susieR)
library(ggplot2)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_prefix <- args[1]
type <- args[2]
N <- as.numeric(args[3])
var_y <- as.numeric(args[4])
coverage <- as.numeric(args[5])
maxit <- as.numeric(args[6])

sdY <- sqrt(var_y)
gwas_file <- paste0(input_prefix, ".locus_results.tsv")
ld_file <- paste0(input_prefix, ".unphased.vcor1")

# Read GWAS dataframe and LD matrix
gwas_df <- fread(gwas_file, sep = "\t")
ld_matrix <- as.matrix(fread(ld_file, sep = "\t", header=FALSE))
rownames(ld_matrix) <- gwas_df$KGP_ID
colnames(ld_matrix) <- gwas_df$KGP_ID

# Prepare coloc data structure
if (var_y > 0) {
    print("Running SuSiE with known outcome variance.")
    coloc_data <- list(
        beta = gwas_df$HARMONIZED_BETA,
        varbeta = gwas_df$SE^2,
        snp = gwas_df$KGP_ID,
        position = gwas_df$POS,
        MAF = gwas_df$MAF,
        N = N,
        type = type,
        sdY = sdY,
        LD = ld_matrix
    )
    results <- runsusie(coloc_data, maxit=maxit, repeat_until_convergence=FALSE, var_y=var_y, n=N, coverage=coverage)
} else {
    print("Running SuSiE with unknown outcome variance.")
    coloc_data <- list(
        beta = gwas_df$HARMONIZED_BETA,
        varbeta = gwas_df$SE^2,
        snp = gwas_df$KGP_ID,
        position = gwas_df$POS,
        MAF = gwas_df$MAF,
        N = N,
        type = type,
        LD = ld_matrix
    )
    results <- runsusie(coloc_data, maxit=maxit, repeat_until_convergence=FALSE, n=N, coverage=coverage)
}

# Make diagnostic plot
eps <- 1e-4
R <- (1 - eps) * coloc_data$LD + eps * diag(nrow(coloc_data$LD))
condz = kriging_rss(coloc_data$beta / sqrt(coloc_data$varbeta), R, n=coloc_data$N)
ggsave(paste0(input_prefix, ".diagnostic_plot.png"), plot = condz$plot, width = 6, height = 4, dpi = 300)

# Make SuSiE PIP plot
png(paste0(input_prefix, ".susie_plot.png"), width = 800, height = 600)
pip_plot <- susie_plot(results, y='PIP', add_legend = TRUE, add_bar = TRUE)
dev.off()

# Make coloc alignement plot
png(paste0(input_prefix, ".coloc_alignement.png"), width = 800, height = 600)
fig <- check_alignment(coloc_data)
dev.off()

# Save SuSiE results
saveRDS(results, file = paste0(input_prefix, ".susie_results.rds"))