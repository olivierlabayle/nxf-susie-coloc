function harmonize_beta(beta, a1, a0, kgp_ref, kgp_alt)
    if (a1 == kgp_alt && a0 == kgp_ref)
        return beta
    elseif (a1 == kgp_ref && a0 == kgp_alt)
        return - beta
    else
        error("Alleles do not match for variant. Cannot harmonize beta.")
    end
end

make_merge_id(c, p, a1, a0) = string(replace(c, "chr" => ""), ":", p, ":", join(sort([a1, a0]), ":"))

string_id_to_merge_id(x) = make_merge_id(split(x, ":")...)

function load_ld_matrix(ld_matrix_file, variants_file)
    return CSV.read(ld_matrix_file, DataFrame; 
        delim='\t', 
        missingstring="NA", 
        header=readlines(variants_file)
    )
end

function compute_LD_matrix(bed_prefix, gwas_variants, locus_output_prefix)
    variants_file = string(locus_output_prefix, ".variants.txt")
    open(variants_file, "w" ) do io
        for var in gwas_variants
            println(io, var)
        end
    end
    run(`plink2 \
        --bfile $bed_prefix \
        --r-unphased square ref-based \
        --extract $(variants_file) \
        --out $locus_output_prefix`
    )
end