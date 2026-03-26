
function write_results_with_merge_ids(results_file, kgp_freqs)
    updated_results_file = "results_and_ref.tsv"
    results = CSV.read(results_file, DataFrame; delim='\t', missingstring="NA")
    select!(results,
        :ID => :GWAS_ID,
        :CHROM,
        :POS,
        :ALLELE_1,
        :ALLELE_0,
        :LOG10P => (x -> abs.(min.(x, 400))) => :LOG10P,
        :BETA,
        :SE,
        :N,
        [:CHROM, :POS, :ALLELE_1, :ALLELE_0] => ByRow(make_merge_id) => :MERGE_ID
    )
    joined_results = innerjoin(results, kgp_freqs, on=:MERGE_ID)
    CSV.write(updated_results_file, joined_results; delim='\t', missingstring="NA")
    return joined_results, updated_results_file
end

function write_significant_clumps(bed_prefix, gwas_results_file;
    min_sig_clump_size = 7,
    lead_pvalue = 5e-8,
    p2_pvalue = 5e-5,
    r2_threshold = 0.2,
    clump_kb = 500,
    clump_id_field = "ID",
    clump_pval_field = "LOG10P",
    allele_1_field = "ALLELE_1"
    )
    output_clump_prefix = "GWAS"
    run(`plink2 --bfile $bed_prefix \
        --clump $gwas_results_file \
        --clump-p1 $lead_pvalue \
        --clump-p2 $p2_pvalue \
        --clump-r2 $r2_threshold \
        --clump-kb $clump_kb \
        --clump-id-field $clump_id_field \
        --clump-log10 \
        --clump-p-field $clump_pval_field \
        --clump-a1-field $allele_1_field \
        --out $output_clump_prefix
    `)
    clumps_file = output_clump_prefix * ".clumps"
    clumps = isfile(clumps_file) ? 
        CSV.read(clumps_file, DataFrame; delim="\t") : 
        DataFrame([col => [] for col in ["#CHROM", "POS", "ID", "NEG_LOG10_P", "TOTAL", "NONSIG", "S0.05", "S0.01", "S0.001", "S0.0001", "SP2"]])
    
    sig_clumps = filter(
        :SP2 => x -> x !== "." && length(split(x, ",")) >= min_sig_clump_size, 
        clumps
    )
    sort!(sig_clumps, ["#CHROM", "POS"])

    CSV.write(
        string(output_clump_prefix, ".sigclumps.tsv"), 
        select(sig_clumps, "#CHROM" => "CHROM", "POS", "ID");
        delim="\t",
        missinstring="NA"
    )

    return 0
end


function prepare_gwas_results(
    results_file,
    kgp_prefix
    )
    # Load KGP frequencies and create merge ids
    @info "Loading Reference frequencies"
    kgp_freqs = CSV.read("$kgp_prefix.afreq", DataFrame; delim='\t', missingstring="NA")
    select!(kgp_freqs,
        "ID" => :KGP_ID,
        "ID" => ByRow(string_id_to_merge_id) => :MERGE_ID,
        "ALT" => :KGP_ALT,
        "REF" => :KGP_REF,
        "ALT_FREQS" => :KGP_ALT_FREQ
    )
    @info "Merging With Reference"
    _, updated_results_file = write_results_with_merge_ids(results_file, kgp_freqs)

    # Get significant clumps
    @info "Finding clumps"
    write_significant_clumps(kgp_prefix, updated_results_file;
        min_sig_clump_size = 7,
        lead_pvalue = 5e-8,
        p2_pvalue = 5e-5,
        r2_threshold = 0.2,
        clump_kb = 500,
        clump_id_field = "KGP_ID",
        clump_pval_field = "LOG10P",
        allele_1_field = "ALLELE_1"
    )
    return 0
end