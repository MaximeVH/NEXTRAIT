### Load modules ###
using DataFrames, CSV, StatsBase, Printf, Plots

### Load data ###
AMP_collection = CSV.read(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\webscraping\recent_amp_scrape\AmP_collection.csv", DataFrame)
AMP_species = CSV.read(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\webscraping\recent_amp_scrape\AmP_species_list.csv", DataFrame)
AMP_Ecoinfo = CSV.read(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\webscraping\recent_amp_scrape\AmP_Ecoinfo.csv", DataFrame)
AMP_tref_parameters = CSV.read(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\webscraping\recent_amp_scrape\AmP_tref_parameters.csv", DataFrame)
AMP_tref_pseudodata = CSV.read(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\webscraping\recent_amp_scrape\AmP_tref_pseudodata.csv", DataFrame)


### Define functions ###
function LifeHistoryTable(AMP_coll)
    unique_vars = unique(AMP_coll.Data)
    AMP_coll_no_guess = AMP_collection[AMP_coll.Reference .!= "guess",:]

    all_counts = []
    all_perc_unique = []
    all_variance = []
    descs = []
    ng_count  = []
    ng_perc_uniqu = []
    ng_varianc = []
    medians = []
    iqrs = []
    q1s = []
    q3s = []
    units_ = []
    means = []

    for (e,i) in enumerate(unique_vars)
        vardf = filter(row -> row.Data == i, AMP_collection)
        var_counts = count(AMP_collection.Data .== i)
        var_perc_unique = length(unique(vardf.Observed))/length(vardf.Observed)
        var_variance = var(vardf.Observed)
        median_ = median(vardf.Observed)
        iqr__ = iqr(vardf.Observed)
        q1 = quantile(vardf.Observed, 0.25)
        q3 = quantile(vardf.Observed, 0.75)
        push!(all_counts, var_counts)
        push!(means, mean(vardf.Observed))
        push!(all_perc_unique, var_perc_unique)
        push!(all_variance, var_variance)
        push!(descs, findmax(countmap(vardf.Description))[2])
        push!(units_, findmax(countmap(vardf.Unit))[2])
        push!(medians, median_)
        push!(iqrs, iqr__)
        push!(q1s, q1)
        push!(q3s, q3)
        vardfng = filter(row -> row.Data == i, AMP_coll_no_guess)
        ng_counts = count(AMP_coll_no_guess.Data .== i)
        ng_perc_unique = length(unique(vardfng.Observed))/length(vardfng.Observed)
        ng_variance = var(vardfng.Observed)
        push!(ng_count, ng_counts)
        push!(ng_perc_uniqu, ng_perc_unique)
        push!(ng_varianc, ng_variance)
    end
    Perc_guess = (1 .- ng_count ./ all_counts) .* 100
    var_perc_unique = all_perc_unique .* 100
    # most_common_element = findmax(countmap(AMP_collection.Description))[2]
    descs = descs .* " (" .* units_ .* ")"
    LifeHistoryDF = DataFrame(Data=unique_vars,Description = descs, Count=all_counts, GuessPerc = Perc_guess, PercUnique=var_perc_unique, Stdev=sqrt.(all_variance),Mean = means, Median=medians, IQR=iqrs)
    sort!(LifeHistoryDF, :Count, rev=true)
    return LifeHistoryDF
end
function print_first_N_rows(df::DataFrame,N=25)
    for row in 1:min(N, nrow(df))
        for col in names(df)
            value = df[row, col]
            if isa(value, Number)
                @printf("%.2f\t", value)
            else
                print(value, "\t")
            end
        end
        println()
    end
end
function DEBParameterTable(AMP_tref_params)
    Params = names(AMP_tref_params[:,9:end])
    Counts = [count(.!isnan.(AMP_tref_params[:,i])) for i in Params]
    Uniques = 100 .* [length(unique(AMP_tref_params[.!isnan.(AMP_tref_params[:,i]),i])) for i in Params] ./ Counts
    Variances = [var(AMP_tref_params[.!isnan.(AMP_tref_params[:,i]),i]) for i in Params]
    stdev = sqrt.(Variances)
    q1s = [quantile(AMP_tref_params[.!isnan.(AMP_tref_params[:,i]),i],0.25) for i in Params]
    q3s = [quantile(AMP_tref_params[.!isnan.(AMP_tref_params[:,i]),i],0.75) for i in Params]
    IQRs = q3s .- q1s
    medians = [median(AMP_tref_params[.!isnan.(AMP_tref_params[:,i]),i]) for i in Params]
    means = [mean(AMP_tref_params[.!isnan.(AMP_tref_params[:,i]),i]) for i in Params]
    DEBParameterDF = DataFrame(Parameter=Params, Count=Counts, Unique=Uniques,Mean=means, Std=stdev, Median=medians, IQR=IQRs)
    sort!(DEBParameterDF, [:Count,:Parameter], rev=[true,false])
    return DEBParameterDF
end
function EstimatedCorrelationHeatmap(AMP_tref_params)
    AMP_tref_params = AMP_tref_params[.!isnan.(AMP_tref_params[:,:E_Hp]),:]
    AMP_tref_params = AMP_tref_params[.!isnan.(AMP_tref_params[:,:E_Hp]),["E_G"  ,"E_Hb" ,"F_m"  ,"T_A"  ,"h_a"  ,"k_J"  ,"kap" ,"kap_P","kap_R","kap_X","p_Am" ,"p_M"  ,"p_T"  ,"s_G"  ,"v" ,"E_Hp"] ]
    sort!(AMP_tref_params, :v, rev=true)
    sorted_col_indices = sortperm(names(AMP_tref_params))
    # Reorder the columns of the DataFrame
    sorted_df = AMP_tref_params[:, sorted_col_indices]
    correlation_matrix = cor(Matrix(sorted_df))
    return heatmap(correlation_matrix, title="Estimated DEB Correlation Matrix", xlabel="Variables", ylabel="Variables", xticks=(1:16, names(sorted_df)), yticks=(1:16, names(sorted_df)),xrotation=45)
end
function MeasuredCorrelationHeatmap(AMP_collect)
    measured_data_considered = ["am","Wwi","Ri","Wwb","Li","ab"]
    unique_amp_spec = unique(AMP_collect.Species)
    counter = 0
    numbs, ams, wwis, ris, wwbs, lis, abs = [], [], [], [], [], [], []
    for spec in unique_amp_spec
        testdf = filter(x-> x.Species == spec, AMP_collect)
        testo = [measured_data_considered[i] in testdf.Data for i in 1:length(measured_data_considered)]
        if all(testo)
            counter += 1
            push!(numbs, sum(testo))
            push!(ams, testdf[testdf.Data .== "am",:Observed][1])
            push!(wwis, testdf[testdf.Data .== "Wwi",:Observed][1])
            push!(ris,testdf[testdf.Data .== "Ri",:Observed][1])
            push!(wwbs, testdf[testdf.Data .== "Wwb",:Observed][1])
            push!(lis, testdf[testdf.Data .== "Li",:Observed][1])
            push!(abs, testdf[testdf.Data .== "ab",:Observed][1])
        end
    end
    EssentialVarsDF = DataFrame(am=ams, Wwi=wwis, Ri=ris, Wwb=wwbs, Li=lis, ab=abs)
    sorted_col_indices = sortperm(names(EssentialVarsDF))
    sorted_df = EssentialVarsDF[:, sorted_col_indices]
    correlation_matrix = cor(Matrix(sorted_df))
    return heatmap(correlation_matrix, title="Measured DEB Correlation Matrix", xlabel="Variables", ylabel="Variables", xticks=(1:6, names(sorted_df)), yticks=(1:6, names(sorted_df)),xrotation=45)
end
function CollectedPerSpecies(AMP_collect)
    measured_data_considered = ["am","Wwi","Ri","Wwb","Li","ab"]
    unique_amp_spec = unique(AMP_collect.Species)
    counter = 0
    specs,numbs, ams, wwis, ris, wwbs, lis, abs = [], [], [], [], [], [], [], []
    for spec in unique_amp_spec
        testdf = filter(x-> x.Species == spec, AMP_collect)
        testo = [measured_data_considered[i] in testdf.Data for i in 1:length(measured_data_considered)]
        if all(testo)
            counter += 1
            push!(numbs, sum(testo))
            push!(specs, spec)
            push!(ams, testdf[testdf.Data .== "am",:Observed][1])
            push!(wwis, testdf[testdf.Data .== "Wwi",:Observed][1])
            push!(ris,testdf[testdf.Data .== "Ri",:Observed][1])
            push!(wwbs, testdf[testdf.Data .== "Wwb",:Observed][1])
            push!(lis, testdf[testdf.Data .== "Li",:Observed][1])
            push!(abs, testdf[testdf.Data .== "ab",:Observed][1])
        end
    end
    EssentialVarsDF = DataFrame(species = specs,am=ams, Wwi=wwis, Ri=ris, Wwb=wwbs, Li=lis, ab=abs)
    return EssentialVarsDF
end
function OverallCorrelationHeatmap(AMP_collect,AMP_tref_params)
    AMP_collect_perspecies = CollectedPerSpecies(AMP_collect)
    IDS = AMP_collect_perspecies.species
    Amp_tref_same_ids = filter(x-> x.ID in IDS, AMP_tref_params)
    params_of_interest = sort(["ID","E_G"  ,"E_Hb" ,"F_m"  ,"T_A"  ,"h_a"  ,"k_J"  ,"kap" ,"kap_P","kap_R","kap_X","p_Am" ,"p_M" ,"s_G"  ,"v" ,"E_Hp"])
    Amp_tref_same_ids_intrest = select(Amp_tref_same_ids, params_of_interest)
    combined_df = hcat(AMP_collect_perspecies, Amp_tref_same_ids_intrest)
    select!(combined_df, Not(:ID))
    sorted_df = select(combined_df, Not(:species))
    # sorted_col_indices = sortperm(names(sorted_df))
    # sorted_df = sorted_df[:, sorted_col_indices]
    correlation_matrix = cor(Matrix(sorted_df))
    return heatmap(correlation_matrix, title="AmP Correlation Matrix", xlabel="Variables", ylabel="Variables", xticks=(1:21, names(sorted_df)), yticks=(1:21, names(sorted_df)),xrotation=45)
end


### Create and display the LifeHistoryTable ###
LifeHistoryDF = LifeHistoryTable(AMP_collection)
print_first_N_rows(LifeHistoryDF)
DEBParameterDF = DEBParameterTable(AMP_tref_parameters)
print_first_N_rows(DEBParameterDF,19)

### Plotting correlation heatmaps ###
EstimatedCorrelationHeatmap(AMP_tref_parameters)
MeasuredCorrelationHeatmap(AMP_collection)
OverallCorrelationHeatmap(AMP_collection,AMP_tref_parameters)
