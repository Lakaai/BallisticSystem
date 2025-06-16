using JSON, DataFrames, CSV, ProgressMeter, ZipFile, Dates
using Base.Filesystem: splitext, mktempdir
using UUIDs
using FilePathsBase

function extract_data(data, df_columns::Dict{String, Any})
    for (k, v) in data
        if isa(v, AbstractDict)  # If the value itself is another dictionary then a key-value pair with raw data has not yet been found
            extract_data(v, df_columns) # Extract the key-values in the nested dictionary
        elseif k == "time_ns" # If the key-value pair is "time_ns" we will handle this case separately 
            df_columns[k] = string(v) # Avoid timestamps being written in scientific notation 
        else 
            df_columns[k] = v # Append the value to the column that has the corresponding key
        end
    end
    return df_columns
end

function export_to_zip(df_dict, filename)
    prefix, _ = splitext(filename)

    zip = ZipFile.Writer("$(prefix)_export.zip") # Create a ZIP file

    for (name, df) in df_dict
        csv_filename = ("$(prefix)_$(name).csv")
        CSV.write(csv_filename, df)
        f = ZipFile.addfile(zip, csv_filename) # Add the CSV file to the ZIP
        write(f, read(csv_filename))
        isfile(csv_filename) && rm(csv_filename)
    end
    close(zip)
end

function readlog_json(filename; data_fields=nothing)

    json_data = JSON.parsefile(filename)

    # If no data_field keyword(s) are provided then return a dictionary that contains a DataFrame for every data-field found in the JSON file
    if data_fields === nothing   

        data_fields = Set{String}() # Container to store the unique data-fields found in the file, these data-fields will each have a unique DataFrame

        # First find all the unique 'data fields' defined by 'polymorphic_name' key 
        for (_, entry) in json_data # Iterate through every data entry in the JSON file
            if haskey(entry, "polymorphic_name") # Check if the entry has the required key
                polymorphic_name = entry["polymorphic_name"] # The 'polymorphic_name' key only appears once so we can assume it is unique
                push!(data_fields, polymorphic_name) # Store the 'polymorphic_name' in data_fields
                println("Data found for struct: ", polymorphic_name)
            end
        end

        # Create the dictionary of DataFrames using the unique data-field strings as keys and DataFrames as values
        df_dict = Dict{String,DataFrame}()

        for (_, entry) in json_data # Iterate over all the JSON data entries
            if  haskey(entry, "ptr_wrapper") && haskey(entry["ptr_wrapper"], "data")
                data = entry["ptr_wrapper"]["data"] # Extract the nested 'data' dictionary

                 # Create a dictionary to store all possible columns for the current 'field' where keys are column names (strings) and values are data of type (any)  
                columns = Dict{String, Any}()

                # Iterate through each 'data' dictionary recursively to extract the key-value pairs and build the columns for the DataFrame
                columns = extract_data(data, columns) 

                # If one of the keys here matches the polymorphic names then create a new data DataFrame for the key
                for (key, _) in data 
                    if key in data_fields && !haskey(df_dict, key)
                        df_dict[key] = DataFrame(columns) # Initialise a DataFrame for the current 'data_field'
                    elseif key in data_fields && haskey(df_dict, key) # If the DataFrame is already constructed then push columns into the existing DataFrame    
                        row = Dict(Symbol(k) => v for (k, v) in columns) # Convert to type Dict{Symbol, Any} because DataFrame expects keys to be type Symbol
                        push!(df_dict[key], row)
                    end 
                end    
            end
        end      
    
        export_to_zip(df_dict, filename) # Export DataFrames to CSV files

    elseif data_fields !== nothing # Handle the case where specific data-field are requested for export
        # Populate the dataframe with the requested data_field(s)
        
        # Create the dictionary of DataFrames using the unique data-field strings as keys and DataFrames as values
        df_dict = Dict{String,DataFrame}()

        for (_, entry) in json_data # Iterate over all the JSON data entries
            if  haskey(entry, "ptr_wrapper") && haskey(entry["ptr_wrapper"], "data")
                data = entry["ptr_wrapper"]["data"] # Extract the nested 'data' dictionary

                 # Create a dictionary to store all possible columns for the current 'field' where keys are column names (strings) and values are data of type (any)  
                columns = Dict{String, Any}()

                # Iterate through each 'data' dictionary recursively to extract the key-value pairs and build the columns for the DataFrame
                columns = extract_data(data, columns) 

                # If one of the keys here matches the polymorphic names then create a new data DataFrame for the key
                for (key, _) in data 
                    if key in data_fields && !haskey(df_dict, key)
                        df_dict[key] = DataFrame(columns) # Initialise a DataFrame for the current 'data_field'
                    elseif key in data_fields && haskey(df_dict, key) # If the DataFrame is already constructed then push columns into the existing DataFrame    
                        row = Dict(Symbol(k) => v for (k, v) in columns) # Convert to type Dict{Symbol, Any} because DataFrame expects keys to be of type Symbol
                        push!(df_dict[key], row)
                    end 
                end    
            end
        end      
        
        export_to_zip(df_dict, filename) # Export DataFrames to CSV files

    end
end

readlog_json("flight_data.json", data_fields=["imu_data_t", "guidance_data_t", "system_data_t"])
# readlog_json("flight_data.json")
