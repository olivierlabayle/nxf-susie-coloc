function cli_settings()
    s = ArgParseSettings(
        description="Susie-Coloc CLI",
        add_version = true,
        commands_are_required = false,
        version=string(pkgversion(PopGen))
    )

    @add_arg_table! s begin
        "finemap-gtex"
            action = :command
            help = "Runs meta-analysis across GWAS results."

    end

    @add_arg_table! s["finemap-gtex"] begin
        "gwas-results"
            arg_type = String
            required = true
            help = "Path to summary statistics."
    
        "--source-software"
            arg_type = String
            default = "saige"
            help = "SOftware that generated the gwas results."

        "--output"
            arg_type = String
            help = "Output filename"
            default = "harmonized_results.tsv"
    end

    return s
end

function julia_main()::Cint
    settings = parse_args(ARGS, cli_settings())
    cmd = settings["%COMMAND%"]
    @info "Running WDL-GWAS CLI: $cmd"
    cmd_settings = settings[cmd]
    if cmd == "make-groups-and-covariates"
        make_groups_and_covariates(
            cmd_settings["covariates-file"];
            groupby_string=cmd_settings["groupby"],
            covariates_string=cmd_settings["covariates"],
            phenotypes_string=cmd_settings["phenotypes"],
            output_prefix=cmd_settings["output-prefix"],
            min_cases_controls=cmd_settings["min-cases-controls"],
            filters_string=cmd_settings["filters"]
        )
    elseif cmd == "merge-covariates-pcs"
        merge_covariates_and_pcs(
            cmd_settings["covariates-file"],
            cmd_settings["pcs-prefix"];
            output=cmd_settings["output"]
        )
    elseif cmd == "make-finemapping-outputs"
        make_finemapping_outputs(
            cmd_settings["merge-list"],
            cmd_settings["gwas-results"];
            output_prefix=cmd_settings["output-prefix"]
        )
    elseif cmd == "make-gwas-outputs"
        make_gwas_outputs(
            cmd_settings["merge-list"];
            output_prefix=cmd_settings["output-prefix"]
        )
    elseif cmd == "finemap"
        finemap_significant_regions(
            cmd_settings["gwas-results-file"],
            cmd_settings["pgen-prefix"],
            cmd_settings["covariates-file"],
            cmd_settings["sample-file"];
            output_prefix=cmd_settings["output-prefix"],
            min_sig_clump_size=cmd_settings["min-sig-clump-size"],
            lead_pvalue=cmd_settings["lead-pvalue"],
            p2_pvalue=cmd_settings["p2-pvalue"],
            r2_threshold=cmd_settings["r2-threshold"],
            clump_kb=cmd_settings["clump-kb"],
            n_causal=cmd_settings["n-causal"],
            phenotype=cmd_settings["phenotype"],
            rss=cmd_settings["rss"],
            susie_max_iter=cmd_settings["susie-max-iter"],
            exclude_string=cmd_settings["exclude"]
        )
    elseif cmd == "meta-analyse"
        meta_analyse(
            cmd_settings["gwas-results-list"];
            output_prefix=cmd_settings["output-prefix"],
            exclude_string=cmd_settings["exclude"],
            method=cmd_settings["method"],
            maf=cmd_settings["maf"]
        )
    elseif cmd == "harmonize-gwas-results"
        harmonize_gwas_results(
            cmd_settings["gwas-results"];
            source_software=cmd_settings["source-software"],
            output=cmd_settings["output"]
        )
    else
        throw(ArgumentError(string("Unknown command: ", cmd)))
    end
    return 0
end