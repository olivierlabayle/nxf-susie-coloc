def get_prefix(file){
    return file.toString().take(file.toString().lastIndexOf('.'))
}

def get_julia_cmd(cpus){
    if (params.USE_SYSIMAGE == false) {
        return "TEMPD=\$(mktemp -d) && JULIA_DEPOT_PATH=\$TEMPD:\$JULIA_DEPOT_PATH julia --project=/opt/FinemapColoc --startup-file=no --threads=${cpus} /opt/FinemapColoc/bin/run.jl"
    }
    else {
        return "TEMPD=\$(mktemp -d) && JULIA_DEPOT_PATH=\$TEMPD:\$JULIA_DEPOT_PATH julia --project=/opt/FinemapColoc --startup-file=no --threads=${cpus} --sysimage=/opt/FinemapColoc/sysimage.so /opt/FinemapColoc/bin/run.jl"
    }        
}

process PrepareGWASResults {
    input:
        path gwas_results
        path ref_files

    output:
        path "results_and_ref.tsv"
        path "GWAS.clumps"

    script:
        def ref_prefix = get_prefix(ref_files[1])
        """
        ${get_julia_cmd(task.cpus)} prepare-gwas-results \
            ${gwas_results} ${ref_prefix}
        """

}

params.GTEX_FILES = "/gpfs/igmmfs01/eddie/ISARIC4C/GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations/*.parquet"
params.GWAS_FINMAPPING_DIRS = "/home/olabayle/isaric/olivier/Covid19/data/gwas_finemapping_results/GWAS_*"
params.GTEX_TISSUES = ["Lung"]
params.GWAS_RESULTS = "/home/olabayle/isaric/olivier/Covid19/data/covid_19_results_2026/meta_analysis_workdir/META_ANALYSIS.all.tsv"
params.REFERENCE_PREFIX = "/gpfs/igmmfs01/eddie/ISARIC4C/olivier/data/kgp-merged-unrelated-or3/kgp.merged.unrelated"

workflow {
    gwas_results = Channel.value(file(params.GWAS_RESULTS))
    reference_files = Channel.fromPath("${params.REFERENCE_PREFIX}*").collect()

    PrepareGWASResults(gwas_results, reference_files)

    // gtex_chr_files_ch = Channel.fromPath(params.GTEX_FILES)
    //     .map { it -> [it.getName().split('\\.'), it] }
    //     .map { it -> [it[0][-2], it[0][0], it[1]] }

    // if (params.GTEX_TISSUES.size() > 0) {
    //     gtex_chr_files_ch = gtex_chr_files_ch.filter { it -> it[1] in params.GTEX_TISSUES }
    // }

    // gwas_chr_lead_pos_ch = Channel.fromPath(params.GWAS_FINMAPPING_DIRS, type: 'dir')
    //     .map { it -> [it.getName().split('_')[1..2], it] }


    // gwas_chr_lead_pos_ch.combine(gtex_chr_files_ch, by: 0).view()
}