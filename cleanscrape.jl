using DataFrames, Gumbo, Cascadia, JSON, HTTP, CSV, Statistics 
url = "https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/species_list.html"
display_table(df) = VSCodeServer.vscodedisplay(df)

function get_species_list_df(url)
    r = HTTP.request("GET", url)
    html_content = String(r.body)
    page = parsehtml(html_content)
    rows = eachmatch(Selector("tr"), page.root)

    df_entries = DataFrame(
    ID = String[],
    Phylum = String[],
    Class = String[],
    Order = String[],
    Family = String[],
    ScientificName = String[],
    CommonName = String[],
    Mod = String[],
    MRE = Float64[],
    SMSE = Float64[],
    COM = Float64[],
    rest = String[],
    )

    for row in rows
        try
            # Extract the ID from the <tr> element
            id = row.attributes["id"]
            
            # Extract the <td> elements within the <tr>
            cells = eachmatch(Selector("td"), row)
            
            # Extract the relevant data from the <td> elements
            phylum = text(cells[1])
            class_ = text(cells[2])
            order = text(cells[3])
            family = text(cells[4])
            scientific_name = text(cells[5])
            common_name = text(cells[6])
            mod = text(cells[7])
            mre = parse(Float64, strip(text(cells[8])))
            smse = parse(Float64, strip(text(cells[9])))
            com = parse(Float64, strip(text(cells[10])))
            d0 = join([text(cells[i]) for i in 11:length(cells)], ", ")
            # d1 = text(cells[18])
            
            # Add the extracted data to the DataFrame
            push!(df_entries, (
                ID = id,
                Phylum = phylum,
                Class = class_,
                Order = order,
                Family = family,
                ScientificName = scientific_name,
                CommonName = common_name,
                Mod = mod,
                MRE = mre,
                SMSE = smse,
                COM = com,
                rest = d0
            ))
        catch e
            # println("Error processing row: ", row)
            continue
        end
    end
    return df_entries
end

function aref_to_texts(aref)
    regex = r"<a\s+href=\"[^\"]*\"\s+title=\"([^\"]*)\">([^<]*)</a>"
    matches = match(regex, aref)
    if matches !== nothing
        title = matches.captures[1]
        value = matches.captures[2]
        return title, value
    else
        return nothing, nothing
    end
end


function get_datatable_from_species(species)
    url = "https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/$(species)/$(species)_res.html"
    r1 = HTTP.request("GET", url)
    html_content = String(r1.body)
    page = parsehtml(html_content)
    links = eachmatch(Selector("table"), page.root)
    infotable = links[2]
    links_colheaders = collect(eachmatch(Selector("th"), infotable))
    links_colheaders = [string(links_colheaders[i][1]) for i in 1:length(links_colheaders)]
    links_data_vals = collect(eachmatch(Selector("td"), infotable))
    links_data_vals = [string(links_data_vals[i][1]) for i in 1:length(links_data_vals)]
    num_cols = length(links_colheaders)
    num_rows = length(links_data_vals) ÷ num_cols
    datamatrix = reshape(links_data_vals,num_cols, num_rows)
    transposed_matrix = permutedims(datamatrix)
    df =  DataFrame(transposed_matrix, links_colheaders)

    newdata = Vector{String}(undef, size(df)[1])
    datanotes = Vector{String}(undef, size(df)[1])
    for (e,i) in enumerate(df.Data)
        if startswith(i, "<a")
            vals = aref_to_texts(i)
            newdata[e] = vals[2]
            datanotes[e] = vals[1]
        else
            newdata[e] = i
            datanotes[e] = ""
        end
    end
    df[!,:Data] = newdata
    df[!,:Datanotes] = datanotes

    newdesc = Vector{String}(undef, size(df)[1])
    descnotes = Vector{String}(undef, size(df)[1])
    for (e,i) in enumerate(df.Description)
        if startswith(i, "<a")
            vals = aref_to_texts(i)
            newdesc[e] = vals[2]
            descnotes[e] = vals[1]
        else
            newdesc[e] = i
            descnotes[e] = ""
        end
    end
    df[!,:Description] = newdesc
    df[!,:Descnotes] = descnotes
    rename!(df,  names(df)[4] => "RE")
    df[!,"RE"] = replace.(replace.(df[:,"RE"], ")" => ""), "(" => "")
    df[!,"Observed"] = parse.(Float64, df[:,:Observed])
    df[!,"Predicted"] = parse.(Float64, df[:,:Predicted])
    df[!,"RE"] = parse.(Float64, df[:,:RE])
    df[!,"Species"] .= species
    return df
end

function get_zero_variate_df(url)
    species_list = get_species_list_df(url);
    IDS = species_list.ID
    df_ = get_datatable_from_species(IDS[1])
    for id in IDS[2:end]
        try
        df_ = vcat(df_, get_datatable_from_species(id))
        catch e
            # println(id)
            continue
        end
    end
    return df_
end

# @time zero_var_df = get_zero_variate_df(url)
# wide_df2 = unstack(zero_var_df, :Species, :Data, :Observed)
# display_table(zero_var_df)

function get_single_pseudodata_df(species,tableind=4)
    url = "https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/$(species)/$(species)_res.html"
    r1 = HTTP.request("GET", url)
    html_content = String(r1.body)
    page = parsehtml(html_content)
    links_pseudodata = eachmatch(Selector("table"), page.root)[tableind]
    links_colheaders = collect(eachmatch(Selector("th"), links_pseudodata))
    links_colheaders = [string(links_colheaders[i][1]) for i in eachindex(links_colheaders)]
    links_data_vals = collect(eachmatch(Selector("td"), links_pseudodata))
    try
        links_data_vals = [string(links_data_vals[i][1]) for i in 1:length(links_data_vals)]
        catch BoundsError
            links_data_vals = [string(links_data_vals[i])[5:end-5] for i in 1:length(links_data_vals)]
    end
    num_cols = length(links_colheaders)
    num_rows = length(links_data_vals) ÷ num_cols
    datamatrix = reshape(links_data_vals,num_cols, num_rows)
    transposed_matrix = permutedims(datamatrix)
    species_name = links_colheaders[3]
    links_colheaders[3] = "Estimated"	
    df_pseudodata =  DataFrame(transposed_matrix, links_colheaders)
    df_pseudodata[!,:Species] .= species_name
    df_pseudodata[!,"Generalised animal"] = parse.(Float64, df_pseudodata[:,"Generalised animal"])
    df_pseudodata[!,"Estimated"] = parse.(Float64, df_pseudodata[:,"Estimated"])
    return df_pseudodata
end

function get_pseudodata_df(url)
    species_list = get_species_list_df(url);
    IDS = species_list.ID
    df_ = get_single_pseudodata_df(IDS[1])
    for id in IDS[2:end]
        try
        df_ = vcat(df_, get_single_pseudodata_df(id))
        catch e
            # println(id)
            continue
        end
    end
    
    return unstack(df_, :Species, :Data, :Estimated)
end

function get_species_params(species) 
    url = "https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/$(species)/$(species)_par.html"
    r1 = HTTP.request("GET", url)
    html_content = String(r1.body)
    page = parsehtml(html_content)
    links = eachmatch(Selector("table"), page.root)
    infotable = links[1]
    links_colheaders = collect(eachmatch(Selector("th"), infotable))
    links_colheaders = [string(links_colheaders[i][1]) for i in 1:length(links_colheaders)]
    links_data_vals = collect(eachmatch(Selector("td"), infotable))
    links_data_vals = [string(links_data_vals[i][1]) for i in 1:length(links_data_vals)]
    links_colheaders = links_colheaders[2:end]
    num_cols = length(links_colheaders)
    num_rows = length(links_data_vals) ÷ num_cols
    datamatrix = reshape(links_data_vals,num_cols, num_rows)
    transposed_matrix = permutedims(datamatrix)
    df =  DataFrame(transposed_matrix, links_colheaders)
    df_new = df[:,1:1]
    df_new[!,"value"] = parse.(Float64, df[:,2])
    return df_new
end

function create_param_dict(df)
    param_dict = Dict()
    for i in 1:size(df)[1]
        param_dict[df[i,1]] = df[i,2]
    end
    return param_dict
end

function get_AMP_parameters_df(url)
    species_list = get_species_list_df(url);
    IDS = species_list.ID
    array_of_dicts = [Dict() for _ in 1:length(IDS)]
    for (n,id) in enumerate(IDS)
        try
            array_of_dicts[n] = create_param_dict(get_species_params(id))
        catch e
            # println(id)
            continue
        end
    end
    keys_ = []
    for dict in array_of_dicts
        for key in keys(dict)
            if !(key in keys_)
                push!(keys_, key)
            end
        end
    end
    param_df = DataFrame()
    for key in keys_
        param_df[!, key] = Float64[]
    end
    parameter_df = species_list[:,1:8]
    for key in keys_
        parameter_df[!, key] = fill(NaN, size(parameter_df)[1])#zeros(Float64, size(parameter_df)[1])
    end
    for (n,dict) in enumerate(array_of_dicts)
        for key in keys_
            if haskey(dict, key)
                parameter_df[n,key] = dict[key]
            end
        end
    end
    return parameter_df
end

function mean_temp(data)
    if all(isnan, data)
        return NaN
    else
        return mean(filter(!isnan, data))
    end
end

function get_AMP_collection_df(url) 
    AMP_spec_list = get_species_list_df(url)
    AMP_collecto = get_zero_variate_df(url)
    datanotes_ = AMP_collecto[:,:Datanotes]
    temps = Vector{Float64}(undef, size(AMP_collecto)[1])
    for (e,note) in enumerate(datanotes_)
        if note == ""
            temps[e] = NaN
        else
            temps[e] = parse(Float64, note[14:end-2])
        end
    end
    AMP_collecto[!,:Temperature] = temps
    AMP_collecto[!, :MeanTemperature] .= NaN
    AMP_collecto[!, :ID] = AMP_collecto.Species
    AMP_collecto2 = leftjoin(AMP_collecto, AMP_spec_list[:,[:ID,:Phylum,:Class,:Order,:Family,:ScientificName]], on=:ID)
    for species in unique(AMP_collecto2.ScientificName)
        # Filter the data for the current species
        species_data = filter(row -> row.ScientificName == species, AMP_collecto2)
        # Calculate the mean temperature for the current species
        mean_temp_value = mean_temp(species_data.Temperature) 
        # Assign the mean temperature to the MeanTemperature column for all entries of the current species
        AMP_collecto2[AMP_collecto2.ScientificName .== species, :MeanTemperature] .= mean_temp_value
    end
    DataNotes_ = []
    DescNotes_ = []
    for i in 1:size(AMP_collecto2)[1]
        if AMP_collecto2.Datanotes[i] == ""
            push!(DataNotes_, missing)
        else
            push!(DataNotes_, AMP_collecto2.Datanotes[i])
        end
        if AMP_collecto2.Descnotes[i] == ""
            push!(DescNotes_, missing)
        else
            push!(DescNotes_, AMP_collecto2.Descnotes[i])
        end
    end
    AMP_collecto2[!,:Datanotes] = DataNotes_
    AMP_collecto2[!,:Descnotes] = DescNotes_
    return AMP_collecto2
end

# Function to extract text from an HTML element
function extract_text(element)
    return isempty(element) ? missing : text(element[1])
end

# Function to extract ecozone information
function extract_ecozone(td)
    ecozone_links = eachmatch(Selector("a"), td)
    ecozones = [text(link) for link in ecozone_links][2:end]
    return join(ecozones, ", ")
end

function parse_html_to_df(ecotable)
    # vals = eachmatch(Selector("td"), ecotable)

    df = DataFrame(
        Model = String[],
        climate = String[],
        migrate = String[],
        phylum = String[],
        COMPLETE = Float64[],
        ecozone = String[],
        food = String[],
        class = String[],
        MRE = Float64[],
        habitat = String[],
        gender = String[],
        order = String[],
        SMSE = Float64[],
        embryo = String[],
        reprod = String[],
        family = String[]
    )

    cells = eachmatch(Selector("td"), ecotable)
    push!(df, (
        # Model = extract_text(eachmatch(Selector("a"), cells[1])),
        Model= extract_ecozone(cells[1]),
        climate = extract_ecozone(cells[2]),
        migrate = extract_ecozone(cells[3]),
        phylum = extract_text(eachmatch(Selector("button"), cells[4])),
        COMPLETE = parse(Float64, strip(split(text(cells[5]), "=")[2])),
        ecozone = extract_ecozone(cells[6]),
        food = extract_ecozone(cells[7]),
        class = extract_text(eachmatch(Selector("button"), cells[8])),
        MRE = parse(Float64, strip(split(text(cells[9]), "=")[2])),
        habitat = extract_ecozone(cells[10]),
        gender = extract_ecozone(cells[11]),
        order = extract_text(eachmatch(Selector("button"), cells[12])),
        SMSE = parse(Float64, strip(split(text(cells[13]), "=")[2])),
        embryo = extract_ecozone(cells[14]),
        reprod = extract_ecozone(cells[15]),
        family = extract_text(eachmatch(Selector("button"), cells[16]))
    ))

    return df
end


function EcoZone_df(species)
    Url_ = "https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/$(species)/$(species)_res.html"
    r1 = HTTP.request("GET", Url_)
    html_content = String(r1.body)
    page = parsehtml(html_content)
    links = eachmatch(Selector("table"), page.root)
    ecotable = links[1]
    df = parse_html_to_df(ecotable)
    df[!,"ID"] .= species
    return df
end

function get_ecozone_df(url)
    species_list = get_species_list_df(url);
    IDS = species_list.ID
    df_ = EcoZone_df(IDS[1])
    for id in IDS[2:end]
        try
        df_ = vcat(df_, EcoZone_df(id))
        catch e
            # println(id)
            continue
        end
    end
    df_reordered = select(df_, :COMPLETE,:MRE,:SMSE, :Model, :climate, :migrate,:ecozone,:habitat, :gender, :food, :embryo, :reprod, :phylum, :class, :order, :family, :ID)
    return df_reordered
end


# @time ecozone_df = get_ecozone_df(url)

# display_table(ecozone_df)