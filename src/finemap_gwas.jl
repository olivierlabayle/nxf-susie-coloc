function finemap_gwas_locus(chrom, lead_pos, kgp_prefix, results; locus_kb = 500, outcome_type="cc", var_y=0.33)
    output_dir = string("gwas_fp_results_chr", chrom, "_", lead_pos)
    isdir(output_dir) || mkdir(output_dir)
    locus_output_prefix = joinpath(output_dir, string("GWAS.chr", chrom, ".", lead_pos))
    @info(string("Fine-mapping locus ", locus_output_prefix))
    # Get locus results
    locus_results = subset(results,
        :CHROM => x -> string.(x) .== string(chrom),
        :POS => x -> x .>= lead_pos - locus_kb * 1000 .&& x .<= lead_pos + locus_kb * 1000,
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

    run(`Rscript $(pkgdir(FinemapColoc))/src/run_susie.R $locus_output_prefix $outcome_type $N $var_y`)

    return 0
end
