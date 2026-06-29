function cli_settings()
    s = ArgParseSettings(
        description="Susie-Coloc CLI",
        add_version = true,
        commands_are_required = false,
        version=string(pkgversion(FinemapColoc))
    )

    @add_arg_table! s begin
        "prepare-gwas-results"
            action = :command
            help = "Prepare GWAS results."

        "finemap-gtex"
            action = :command
            help = "Finemap GTEX results."
        
        "finemap-gwas-locus"
            action = :command
            help = "Finemap GWAS locus."

        "aggregate-coloc-results"
            action = :command
            help = "Aggregate coloc results."
    end

    @add_arg_table! s["prepare-gwas-results"] begin
        "results-file"
            arg_type = String
            help = "Path to GWAS summary statistics"

        "ref-prefix"
            arg_type = String
            help = "Prefix to reference dataset"

        "--min-sig-clump-size"
            arg_type = Int
            help     = "Minimum number of SNPs in a clump to be considered significant"
            default  = 7
        
        "--lead-pvalue"
            arg_type = Float64
            help     = "Lead p-value threshold"
            default  = 5e-8

        "--p2-pvalue"
            arg_type = Float64
            help     = "Secondary p-value threshold"
            default  = 5e-5

        "--r2-threshold"
            arg_type = Float64
            help     = "R2 threshold for clumping"
            default  = 0.2
        
        "--clump-kb"
            arg_type = Int
            help     = "Clumping distance in kb"
            default  = 500
    end

    @add_arg_table! s["finemap-gtex"] begin
        "gtex-file"
            arg_type = String
            help = "Path to parquet GTEX file."
    
        "gwas-locus-results-dir"
            arg_type = String
            help = "Path to locus GWAS finemapping results directory."

        "chrom"
            arg_type = String
            help = "Chromosome"
        
        "lead-pos"
            arg_type = Int
            help     = "Lead SNP position in the locus"
        
        "tissue"
            arg_type = String
            help = "GTEX Tissue"

        "N"
            arg_type = Int
            help     = "Number of samples used by the GTEX analysis"
        
        "--coverage"
            arg_type = Float64
            help     = "SuSiE coverage parameter"
            default  = 0.95

        "--susie-maxit"
            arg_type = Int
            help     = "Maximum SuSiE iterations."
            default = 1000
    end

    @add_arg_table! s["finemap-gwas-locus"] begin
        "chrom"
            arg_type = String
            help = "Chromosome of the lead variant."

        "pos"
            arg_type = Int
            help = "Position of the lead variant."

        "ref-prefix"
            arg_type = String
            help = "Prefix to reference dataset."

        "gwas-results"
            arg_type = String
            help = "Path to prepare gwas results."

        "--locus-kb"
            arg_type = Int
            help     = "Lead SNP position in the locus"
            default = 500
        
        "--outcome-type"
            arg_type = String
            help = "Type of GWAS outcome cc/quant."
            default = "cc"
        
        "--var-y"
            arg_type = Float64
            help     = "Variance of the outcome if known"
            default  = 0.0

        "--coverage"
            arg_type = Float64
            help     = "SuSiE coverage parameter"
            default  = 0.95

        "--susie-maxit"
            arg_type = Int
            help     = "Maximum SuSiE iterations."
            default = 1000
    end

    @add_arg_table! s["aggregate-coloc-results"] begin
        "results-file-list"
            arg_type = String
            help = "Path to list of files to aggregate"

        "--output"
            arg_type = String
            default = "aggregated_coloc_results.tsv"
            help = "Output path"
    end

    return s
end

function julia_main()::Cint
    settings = parse_args(ARGS, cli_settings())
    cmd = settings["%COMMAND%"]
    @info "Running SuSiE-Coloc CLI: $cmd"
    cmd_settings = settings[cmd]
    if cmd == "finemap-gtex"
        finemap_gtex_file(
            cmd_settings["gtex-file"],
            cmd_settings["gwas-locus-results-dir"],
            cmd_settings["chrom"],
            cmd_settings["lead-pos"],
            cmd_settings["tissue"],
            cmd_settings["N"];
            coverage=cmd_settings["coverage"],
            susie_maxit=cmd_settings["susie-maxit"]
        )
    elseif cmd == "prepare-gwas-results"
        prepare_gwas_results(
            cmd_settings["results-file"],
            cmd_settings["ref-prefix"];
            min_sig_clump_size = cmd_settings["min-sig-clump-size"],
            lead_pvalue = cmd_settings["lead-pvalue"],
            p2_pvalue = cmd_settings["p2-pvalue"],
            r2_threshold = cmd_settings["r2-threshold"],
            clump_kb = cmd_settings["clump-kb"]
        )
    elseif cmd == "finemap-gwas-locus"
        finemap_gwas_locus(
            cmd_settings["chrom"],
            cmd_settings["pos"],
            cmd_settings["ref-prefix"],
            cmd_settings["gwas-results"];
            locus_kb = cmd_settings["locus-kb"],
            outcome_type = cmd_settings["outcome-type"],
            var_y=cmd_settings["var-y"],
            coverage=cmd_settings["coverage"],
            susie_maxit=cmd_settings["susie-maxit"]
        )
    elseif cmd == "aggregate-coloc-results"
        aggregate_coloc_results(
            cmd_settings["results-file-list"];
            output=cmd_settings["output"]
        )
    else
        throw(ArgumentError(string("Unknown command: ", cmd)))
    end
    return 0
end