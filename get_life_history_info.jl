include("cleanscrape.jl")

@time begin
AmP_species_list = get_species_list_df(url)
CSV.write("AmP_species_list.csv", AmP_species_list)
AmP_pseudodata = get_pseudodata_df(url)
CSV.write("AmP_tref_pseudodata.csv", AmP_pseudodata)
AmP_tref_parameters = get_AMP_parameters_df(url)
CSV.write("AmP_tref_parameters.csv", AmP_tref_parameters)
AmP_collection = get_AMP_collection_df(url) 
CSV.write("AmP_collection.csv", AmP_collection)
AmP_EcoInfo = get_ecozone_df(url)
CSV.write("AmP_EcoInfo.csv", AmP_EcoInfo)
end

