function finemap_gwas_locus(clump, kgp_prefix, results, output_dir; locus_kb = 500, var_y=0.33)
    locus_dir = joinpath(output_dir, string("GWAS_chr", clump["#CHROM"], "_", clump["POS"]))
    isdir(locus_dir) || mkdir(locus_dir)
    locus_output_prefix = joinpath(locus_dir, string("GWAS.chr", clump["#CHROM"], ".", clump["POS"]))
    if isfile(string(locus_output_prefix, ".susie_results.rds"))
        @info(string("Fine-mapping results for locus ", locus_output_prefix, " already exist. Skipping fine-mapping."))
        return
    end
    @info(string("Fine-mapping locus ", locus_output_prefix))
    # Get locus results
    locus_results = subset(results,
        :CHROM => x -> string.(x) .== string(clump["#CHROM"]),
        :POS => x -> x .>= clump["POS"] - locus_kb * 1000 .&& x .<= clump["POS"] + locus_kb * 1000,
        :KGP_ALT_FREQ => x -> x .!= 0 .&& x .!= 1
    )
    # Harmonize the effect alleles and betas with the KGP reference
    transform!(locus_results, 
        [:BETA, :ALLELE_1, :ALLELE_0, :KGP_REF, :KGP_ALT] => ByRow(harmonize_beta) => :HARMONIZED_BETA,
        :KGP_ALT_FREQ => ByRow(x -> x >= 0.5 ? 1 .- x : x) => :MAF
    )
    # Sort locus results by position and write to file
    sort!(locus_results, :POS)
    CSV.write(
        string(locus_output_prefix, ".locus_results.tsv"), 
        locus_results; 
        delim='\t', 
        missingstring="NA"
    )
    # Compute LD matrix for locus variants an check alignement
    compute_LD_matrix(kgp_prefix, locus_results.KGP_ID, locus_output_prefix)
    readlines(string(locus_output_prefix, ".unphased.vcor1.vars")) == locus_results.KGP_ID || error("Variants in LD matrix do not match locus results.")

    N = Int(median(locus_results.N))

    run(`Rscript src/susie_finemap.R $locus_output_prefix $N $var_y`)

    return locus_results
end
