function finemap_gtex_locus(gene_qtl_df, ld_matrix, gene_output_dir; N=940)
    isdir(gene_output_dir) || mkdir(gene_output_dir)
    locus_output_prefix = joinpath(gene_output_dir, gene_output_dir)

    @info(string("Fine-mapping locus ", gene_output_dir))
    # Harmonize the effect alleles and betas with the KGP reference
    transform!(gene_qtl_df, 
        [:BETA, :ALLELE_1, :ALLELE_0, :KGP_REF, :KGP_ALT] => ByRow(harmonize_beta) => :HARMONIZED_BETA,
    )
    # Sort locus results by position and write to file
    sort!(gene_qtl_df, :POS)
    CSV.write(
        string(locus_output_prefix, ".locus_results.tsv"), 
        gene_qtl_df; 
        delim='\t',
        missingstring="NA"
    )
    # Compute LD matrix for locus variants an check alignement
    variants_to_keep_indices = indexin(gene_qtl_df.KGP_ID, names(ld_matrix))
    ld_matrix_to_keep = ld_matrix[variants_to_keep_indices, variants_to_keep_indices]

    CSV.write(string(locus_output_prefix, ".unphased.vcor1"), ld_matrix_to_keep; delim='\t', header=false)

    run(`Rscript src/susie_finemap.R $locus_output_prefix quant $N 0`)

    return 0
end

"""
Variants are encoded as: chr_pos_ref_alt_assembly.
The alt allele is the effect allele in GTEx, which is annotated as ALLELE_1 while the ref allele is annotated as ALLELE_0.
"""
function extract_gtex_info_from_id(variant_id)
    chr, pos, a0, a1, assembly = split(variant_id, "_")
    return replace(chr, "chr" => ""), parse(Int, pos), a1, a0
end

function finemap_gtex_file(
    gtex_file,
    gwas_locus_results_dir,
    chrom,
    lead_pos,
    tissue,
    N
    )
    
    # Find relevant files
    gwas_dir_files = readdir(gwas_locus_results_dir, join=true)
    gwas_locus_results_file = gwas_dir_files[findfirst(endswith(".locus_results.tsv"), gwas_dir_files)]
    ld_matrix_file = gwas_dir_files[findfirst(endswith(".unphased.vcor1"), gwas_dir_files)]
    variants_file = gwas_dir_files[findfirst(endswith(".unphased.vcor1.vars"), gwas_dir_files)]

    # Load GWAS fine-mapping results for the locus of interest
    @info "Loading GWAS fine-mapping results for locus $chrom:$lead_pos"
    gwas_locus_results = CSV.read(
        gwas_locus_results_file, 
        DataFrame; 
        delim='\t', 
        missingstring="NA",
        select=[:MERGE_ID, :KGP_ID, :MAF, :KGP_ALT, :KGP_REF]
    )
    # Load GTEx results for the tissue and chromosome of interest
    @info "Loading GTEx Data for tissue $tissue and chromosome $chrom"
    gtex_ds = Parquet2.Dataset(gtex_file)
    gtex_df = DataFrame(gtex_ds; copycols=false)
    transform!(gtex_df,
        :variant_id => ByRow(extract_gtex_info_from_id) => [:CHROM, :POS, :ALLELE_1, :ALLELE_0],
    )
    transform!(gtex_df,
        [:CHROM, :POS, :ALLELE_1, :ALLELE_0] => ByRow(make_merge_id) => :MERGE_ID
    )
    joined_results = innerjoin(gtex_df, gwas_locus_results, on=:MERGE_ID)
    rename!(joined_results, :slope => :BETA, :slope_se => :SE, :pval_nominal => :PVALUE)

    # Load LD matrix for the locus variants
    ld_matrix = load_ld_matrix(ld_matrix_file, variants_file)

    # Finemap each gene expression for the locus
    for (gene_key, gene_qtl_df) in pairs(groupby(joined_results, :gene_id))
        @info string("Fine-mapping gene ", gene_key.gene_id)
        gene_output_dir = string("GTEX_", chrom, "_", lead_pos, "_", tissue, "_", gene_key.gene_id)
        finemap_gtex_locus(gene_qtl_df, ld_matrix, gene_output_dir; N=N)
    end

    return 0
end