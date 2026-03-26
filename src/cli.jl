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
    end

    @add_arg_table! s["prepare-gwas-results"] begin
        "results-file"
            arg_type = String
            help = "Path to GWAS summary statistics"

        "ref-prefix"
            arg_type = String
            help = "Prefix to reference dataset"

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
            cmd_settings["N"]
        )
    elseif cmd == "prepare-gwas-results"
        prepare_gwas_results(
            cmd_settings["results-file"],
            cmd_settings["ref-prefix"]
        )
    elseif cmd == "finemap-gwas-locus"
        finemap_gwas_locus(
            cmd_settings["chrom"],
            cmd_settings["pos"],
            cmd_settings["ref-prefix"],
            cmd_settings["gwas-results"];
            locus_kb = cmd_settings["locus-kb"],
            outcome_type = cmd_settings["outcome-type"],
            var_y=cmd_settings["var-y"]
        )
    else
        throw(ArgumentError(string("Unknown command: ", cmd)))
    end
    return 0
end