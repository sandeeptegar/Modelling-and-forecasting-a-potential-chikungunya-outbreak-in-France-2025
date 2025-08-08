function read_case_files(start_date::Date, end_date::Date)
    all_data = DataFrame[]  # Vector to collect DataFrames

    for date in start_date:end_date
        # Construct the expected date part in the filename (e.g., "dc_21June")
        date_prefix = "dc_$(Dates.format(date, "dU"))"

        # Look inside the simout directory
        for file in readdir("simout")
            if occursin(date_prefix, file) && endswith(file, ".csv")
                full_path = joinpath("simout", file)
                df = CSV.read(full_path, DataFrame)
                push!(all_data, df)
            end
        end
    end

    return all_data
end
