params.GTEX_FILES = "/gpfs/igmmfs01/eddie/ISARIC4C/GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations/*.parquet"
params.GWAS_FINMAPPING_DIRS = "/home/olabayle/isaric/olivier/Covid19/data/gwas_finemapping_results/GWAS_*"
params.GTEX_TISSUES = ["Lung"]

workflow {
    gtex_chr_files_ch = Channel.fromPath(params.GTEX_FILES)
        .map { it -> [it.getName().split('\\.'), it] }
        .map { it -> [it[0][-2], it[0][0], it[1]] }

    if (params.GTEX_TISSUES.size() > 0) {
        gtex_chr_files_ch = gtex_chr_files_ch.filter { it -> it[1] in params.GTEX_TISSUES }
    }

    gwas_chr_lead_pos_ch = Channel.fromPath(params.GWAS_FINMAPPING_DIRS, type: 'dir')
        .map { it -> [it.getName().split('_')[1..2], it] }


    gwas_chr_lead_pos_ch.combine(gtex_chr_files_ch, by: 0).view()
}