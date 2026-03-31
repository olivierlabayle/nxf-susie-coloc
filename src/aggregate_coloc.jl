function aggregate_coloc_results(results_file_list; output="aggregated_coloc_results.tsv")
    result_files = readlines(results_file_list)
    results = mapreduce(vcat, result_files) do filename
        CSV.read(filename, DataFrame)
    end
    CSV.write(output, results; delim="\t")
end