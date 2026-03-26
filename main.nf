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
    label 'multithreaded'
    label 'mediummem'
    publishDir 'results/prepare_gwas_results'

    input:
        path gwas_results
        path ref_files

    output:
        path "results_and_ref.tsv", emit: gwas_results
        path "GWAS.clumps", emit: clumps

    script:
        def ref_prefix = get_prefix(ref_files[1])
        """
        ${get_julia_cmd(task.cpus)} prepare-gwas-results \
            ${gwas_results} ${ref_prefix}
        """
}

process FinemapGWASLocus {
    label 'multithreaded'
    label 'mediummem'
    publishDir 'results/gwas_fp_results'

    input:
        tuple val(chrom), val(pos)
        path(ref_files)
        path(gwas_results)

    output:
        tuple val(chrom), val(pos), path("gwas_fp_results_chr${chrom}_${pos}", type: 'dir')

    script:
    def ref_prefix = get_prefix(ref_files[1])
    """
    ${get_julia_cmd(task.cpus)} finemap-gwas-locus \
        ${chrom} \
        ${pos} \
        ${ref_prefix} \
        ${gwas_results} \
        --locus-kb=${params.LOCUS_KB} \
        --outcome-type=${params.GWAS_OUTCOME_TYPE} \
        --var-y=${params.GWAS_OUTCOME_VAR}
    """
}

workflow {
    gwas_results = Channel.value(file(params.GWAS_RESULTS))
    reference_files = Channel.fromPath("${params.REFERENCE_PREFIX}*").collect()

    prepared_results = PrepareGWASResults(gwas_results, reference_files)

    clumps_ch = prepared_results.clumps
        .splitCsv(sep: "\t", skip: 1, header: ["CHROM", "POS", "ID", "NEG_LOG10_P", "TOTAL", "NONSIG", "S0.05",	"S0.01", "S0.001", "S0.0001", "SP2"])
        .map { it -> [it.CHROM, it.POS]}

    gwas_fp_ch = FinemapGWASLocus(clumps_ch, reference_files, prepared_results.gwas_results)

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