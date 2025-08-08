function seeding_run(seeding_start::Date, seeding_end::Date)
    # Ensure seeding_start is before or equal to seeding_end
    if seeding_start > seeding_end
        throw(ArgumentError("seeding_start should be less than or equal to seeding_end"))
    end
    
    # Loop from seeding_start to seeding_end inclusive
    for date in seeding_start:seeding_end
        # simulate model
        # Define start date as a string
        start_date_str = "18.11.2021"
        seed_date_str = string(Dates.format(Date(date), "dd.mm.yyyy"))
        sol1 = chikungunya(43.77, 4.29, start_date_str, "30.12.2025", seed_date_str)
        
        # Get daily cases
        daily_case = diff(sol1["i_inst"] .+ sol1["Hr"])
        # Trim Hr to match diff length (diff returns n-1 elements)
        hr_trimmed = sol1["Hr"][2:end]

        # Create DataFrame with both daily_case and corresponding Hr
        output = DataFrame(
            time        = DataFrame(daily_case).timestamp,
            daily_cases = DataFrame(daily_case).i_inst_Hr,
            cum_cases   = DataFrame(hr_trimmed).Hr
            )
        cumsum = Int(round(maximum(values(sol1.Hr))))
        
        # Write output to file with interpolated date
        # Generate the filename
        filename = "simout/dc_$(Dates.format(date, "dU"))_$(cumsum).csv"
        CSV.write(filename, output)
    end
end