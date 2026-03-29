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
        path "GWAS.sigclumps.tsv", emit: sig_clumps

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
        tuple val(chrom), val(pos), path("gwas_fp_results_chr${chrom}_${pos}", type: 'dir'), eval("cat gwas_fp_results_chr${chrom}_${pos}/status.txt")

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
        --var-y=${params.GWAS_OUTCOME_VAR} \
        --coverage=${params.SUSIE_COVERAGE} \
        --susie-maxit=${params.SUSIE_MAXIT}
    """
}

process FinemapGTEXTFile {
    label 'multithreaded'
    label 'mediummem'
    publishDir 'results/gtex_fp_results'

    input:
        tuple val(chrom), val(pos), path(gwas_fp_dir), val(tissue), path(gtex_file)

    output:
        tuple val(chrom), val(pos), path("gwas_fp_results_chr${chrom}_${pos}", type: 'dir')

    script:
    def ref_prefix = get_prefix(ref_files[1])
    """
    ${get_julia_cmd(task.cpus)} finemap-gtex \
        ${gtex_file} \
        ${gwas_fp_dir} \
        ${chrom} \
        ${pos} \
        ${tissue} \
        ${params.GTEX_SAMPLE_SIZE} \
        --coverage=${params.SUSIE_COVERAGE} \
        --susie-maxit=${params.SUSIE_MAXIT}
    """
}

workflow {
    gwas_results = Channel.value(file(params.GWAS_RESULTS))
    reference_files = Channel.fromPath("${params.REFERENCE_PREFIX}*").collect()

    // Prepare GAWS results by Merging with the reference datasets and identifying significant loci
    prepared_results = PrepareGWASResults(gwas_results, reference_files)

    // Finemap significant GWAS loci
    sig_clumps_ch = prepared_results.sig_clumps
        .splitCsv(sep: "\t", skip: 1, header: ["CHROM", "POS", "ID"])
        .map { it -> [it.CHROM, it.POS]}
    gwas_fp_ch = FinemapGWASLocus(sig_clumps_ch, reference_files, prepared_results.gwas_results)

    gwas_fp_ch.view()
    // // Read GTEX files
    // gtex_chr_files_ch = Channel.fromPath(params.GTEX_FILES)
    //     .map { it -> [it.getName().split('\\.'), it] }
    //     .map { it -> [it[0][-2], it[0][0], it[1]] }
    // if (params.GTEX_TISSUES.size() > 0) {
    //     gtex_chr_files_ch = gtex_chr_files_ch.filter { it -> it[1] in params.GTEX_TISSUES }
    // }

    // // Finemap GWAS matching loci GTEX results
    // FinemapGTEXTFile(gwas_fp_ch.combine(gtex_chr_files_ch, by: 0))
}