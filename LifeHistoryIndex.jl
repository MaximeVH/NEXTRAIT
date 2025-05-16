using CSV, DataFrames, Statistics
using MLJ, Plots
using MLJMultivariateStatsInterface
using MultivariateStats
AMPCOL = CSV.read(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\webscraping\recent_amp_scrape\AmP_collection.csv", DataFrame)
display_table(df) = VSCodeServer.vscodedisplay(df)
function convert_columns_to_float64(df::DataFrame)
    for col in names(df)[2:end]
        df[!, col] = convert(Vector{Float64}, df[!, col])
    end
    return df
end
function orders_of_magnitude_covered(values::Vector{Float64})
    orders_of_magnitude = log10.(values)
    min_order = floor(minimum(orders_of_magnitude))
    max_order = floor(maximum(orders_of_magnitude))
    num_orders_of_magnitude = max_order - min_order
    return num_orders_of_magnitude
end
function standardize_columns(df::DataFrame)
    for col in names(df)[2:end]
        col_mean = mean(df[!, col])
        col_std = std(df[!, col])
        df[!, col] = (df[!, col] .- col_mean) ./ col_std
    end
        return df
end
function standardize_columns_multi(dftrain, dftest)
    for col in names(dftrain)[2:end]
        col_mean = mean(dftrain[:, col])
        col_std = std(dftrain[:, col])
        dftrain[:, col] = (dftrain[:, col] .- col_mean) ./ col_std
        dftest[:, col] = (dftest[:, col] .- col_mean) ./ col_std
    end
    return dftrain, dftest
end
function GetLifeHistoryDFOld(AMP_collect,removeguess=false) #can be adjusted to include a set of chosen variables, where a DF is created with as many columns as variables 
    measured_data_considered = ["am","Wwi","Ri","Wwb","Li","ab"]
    if removeguess
        AMP_collect = filter(x -> x.Reference != "guess", AMP_collect)
    end
    unique_amp_spec = unique(AMP_collect.Species)
    counter = 0
    numbs, ams, wwis, ris, wwbs, lis, abs = [], [], [], [], [], [], []
    species = []
    for spec in unique_amp_spec
        testdf = filter(x-> x.Species == spec, AMP_collect)
        testo = [measured_data_considered[i] in testdf.Data for i in 1:length(measured_data_considered)]
        if all(testo)
            counter += 1
            push!(numbs, sum(testo))
            push!(species, spec)
            push!(ams, testdf[testdf.Data .== "am",:Observed][1])
            push!(wwis, testdf[testdf.Data .== "Wwi",:Observed][1])
            push!(ris,testdf[testdf.Data .== "Ri",:Observed][1])
            push!(wwbs, testdf[testdf.Data .== "Wwb",:Observed][1])
            push!(lis, testdf[testdf.Data .== "Li",:Observed][1])
            push!(abs, testdf[testdf.Data .== "ab",:Observed][1])
        end
    end
    # EssentialVarsDF = DataFrame(species = species, am=ams, Wwi=wwis, Ri=ris, Wwb=wwbs, Li=lis, ab=abs)
    EssentialVarsDF = DataFrame(species = species, am=log10.(ams), Wwi=log10.(wwis), Ri=log10.(ris), Wwb=log10.(wwbs), Li=log10.(lis), ab=log10.(abs))
    return standardize_columns(convert_columns_to_float64(EssentialVarsDF))
end
function GetLifeHistoryDF(AMP_collect;vars=["am","Wwi","Ri","Wwb","Li","ab"],removeguess=false) 
    if removeguess
        AMP_collect = filter(x -> x.Reference != "guess", AMP_collect)
    end
    unique_amp_spec = unique(AMP_collect.Species)
    varsinfo = [[] for i in 1:length(vars)]
    species = []
    for spec in unique_amp_spec
        testdf = filter(x-> x.Species == spec, AMP_collect)
        testo = [vars[i] in testdf.Data for i in 1:length(vars)]
        if all(testo)
            push!(species, spec)
            for (e,v) in enumerate(vars)
                push!(varsinfo[e], testdf[testdf.Data .== v,:Observed][1])
            end
        end
    end
    EssentialVarsDF = hcat(DataFrame(species = species),DataFrame(log10.(hcat(varsinfo...)), vars))
    return standardize_columns(convert_columns_to_float64(EssentialVarsDF))
end
function LifeHistoryPCA(LH_df,var="Wwi")
    X = select(LH_df, Not(:species))
    # Perform PCA
    PCA = @load PCA pkg=MultivariateStats
    pca_model = PCA(maxoutdim=2)
    mach = machine(pca_model, X)
    fit!(mach)
    # Transform the data using the PCA model
    X_pca = MLJ.transform(mach, X)
    explained_vars = explained_var(Matrix(X))
    # Plot the first two components
    pcafig = scatter(X_pca[:, 1], X_pca[:, 2], marker_z=LH_df[:,var],color=:viridis ,xlabel="PC1 ($(explained_vars[1])%)", ylabel="PC2 ($(explained_vars[2])%)", title="Life History PCA, colored by $var", colorbar=true,markerstrokewidth=0)
    return pcafig,X_pca
end
function MitoGenomeLHIRegression(mito_df,LH_df;kfolds=6,kmers=4,model=LightGBM.MLJInterface.LGBMRegressor)
    PCA = @load PCA pkg=MultivariateStats
    pipe_regression = ContinuousEncoder |> Standardizer |> model
    # step one: get the mitogenomes that have life history info
    mito_df = filter(x -> x.ID in LH_df.species, mito_df)
    # step two: get the kmer features
    seqfeats = generate_kmer_features(remove_non_ATGC.(mito_df.sequence),kmers)
    # step three: get the folds
    Folds = species_stratified_folds(mito_df.ID, kfolds)
    y_tests = []
    y_preds = []
    for (traininds,testinds) in Folds #
        # step four: standardisation and PCA of training data
        dftrain_Z, dftest_Z = standardize_columns_multi(LH_df[traininds,:], LH_df[testinds,:])
        Ytrain = select(dftrain_Z, Not(:species))
        Ytest = select(dftest_Z, Not(:species))
        # Transform the data using the PCA model
        pca_model = PCA(maxoutdim=2)
        machPCA = machine(pca_model, Ytrain)
        fit!(machPCA)
        Ytrain_pca = MLJ.transform(machPCA, Ytrain)
        Ytest_pca = MLJ.transform(machPCA, Ytest)
        ytrain = Ytrain_pca[:,1]
        ytest = Ytest_pca[:,1]
        Xtrain = seqfeats[traininds,:]
        Xtest = seqfeats[testinds,:]
        RegModel = model()
        machReg = machine(pipe_regression, Xtrain, ytrain)
        MLJ.fit!(machReg, verbosity=0)
        ypred = MLJ.predict(machReg, Xtest)
        y_tests = vcat(y_tests, ytest)
        y_preds = vcat(y_preds, ypred)
    end

    # step six: regression and evaluation
    # X,y = dropnans(seqfeats, vardf.Observed)
    # RegModel = model()
    # pipe_regression = ContinuousEncoder |> Standardizer |> RegModel
    # kmermito_results = evaluate(pipe_regression, X, y, resampling=CV(shuffle=true,nfolds=6),
    # measures=[l1,l2,rms,rsq_], verbosity=0)
    # return kmermito_results
    return y_tests, y_preds
end
function rsquared(y_true::Vector{Float64}, y_pred::Vector{Float64})
    ss_res = sum((y_true .- y_pred).^2)
    ss_tot = sum((y_true .- mean(y_true)).^2)
    r2 = 1 - ss_res / ss_tot
    return r2
end
function GetHabitatTable(EcoInfo)
    life_stage_dict = Dict('0' => 1,'o'=>1,'O'=>1, 'b' => 2, 'x' => 3, 'j' => 4, 'p' => 5, 'e' => 6, 'i' => 7)
    LifeStageTable = fill("", size(EcoInfo,1), 7)
    habitatinfo = EcoInfo.habitat
    for i in 1:length(habitatinfo)
        hinfo = split(replace(habitatinfo[i], " " => ""), ",")
        if length(hinfo) == 1
            firstind = life_stage_dict[hinfo[1][1]]
            lastind = life_stage_dict[hinfo[1][2]] - 1
            LifeStageTable[i,firstind:lastind] .= hinfo[1][3:end]
        else
            for j in 1:length(hinfo)
                firstind = life_stage_dict[hinfo[j][1]]
                lastind = life_stage_dict[hinfo[j][2]] - 1
                LifeStageTable[i,firstind:lastind] .= hinfo[j][3:end]
            end
        end
    end
    LifeStageTable = DataFrame(LifeStageTable, ["startdev", "birth", "weaning", "juvenile", "adult", "pupa", "death"])
    prevnames = names(LifeStageTable)
    LifeStageTable[!,:ID] = EcoInfo.ID
    return select(LifeStageTable, ["ID", prevnames...])
end
function GetFoodTable(EcoInfo)
    life_stage_dict = Dict('0' => 1,'o'=>1,'O'=>1, 'b' => 2, 'x' => 3, 'j' => 4, 'p' => 5, 'e' => 6, 'i' => 7)
    LifeStageTable = fill("", size(EcoInfo,1), 7)
    foodinfo = EcoInfo.food
    for i in 1:length(foodinfo)
        finfo = split(replace(foodinfo[i], " " => ""), ",")
        if length(finfo) == 1
            firstind = life_stage_dict[finfo[1][1]]
            lastind = life_stage_dict[finfo[1][2]]-1
            LifeStageTable[i,firstind:lastind] .*= finfo[1][3:end]
        else
            for j in 1:length(finfo)
                firstind = life_stage_dict[finfo[j][1]]
                lastind = life_stage_dict[finfo[j][2]]-1
                LifeStageTable[i,firstind:lastind] .*= finfo[j][3:end]
            end
        end
    end
    LifeStageTable = DataFrame(LifeStageTable, ["startdev", "birth", "weaning", "juvenile", "adult", "pupa", "death"])
    prevnames = names(LifeStageTable)
    LifeStageTable[!,:ID] = EcoInfo.ID
    return select(LifeStageTable, ["ID", prevnames...])
end
# rego = NearestNeighborModels.KNNRegressor
function simplifyLifeStageTable!(table)
    for col in names(table)[2:end]
        for row in 1:size(table,1)
            table[row,col] = join(unique(replace(table[row,col], r"[^A-Z]" => "")))
        end
    end
    return table 
end
function SimpleEcozoneInfo(EcoInfo)
    simple_ecozone = [string(EcoInfo.ecozone[i][1]) for i in 1:size(EcoInfo,1)]
    simple_gender = [string(EcoInfo.gender[i][1]) for i in 1:size(EcoInfo,1)]
    simple_reproduction = [string(EcoInfo.reprod[i][1]) for i in 1:size(EcoInfo,1)]
    simple_embryo = [string(EcoInfo.embryo[i][1]) for i in 1:size(EcoInfo,1)]
    BasicEcoInfo = DataFrame(ID = EcoInfo.ID, ecozone = simple_ecozone, gender = simple_gender, reproduction = simple_reproduction, embryo = simple_embryo)
    food_table = GetFoodTable(EcoInfo)
    simplifyLifeStageTable!(food_table)
    habitat_table = GetHabitatTable(EcoInfo)
    simplifyLifeStageTable!(habitat_table)
    return BasicEcoInfo, food_table, habitat_table
end
function explained_var(X,outdim=2)
    M = fit(MultivariateStats.PCA, X, maxoutdim=outdim)
    var_explained = Int.(round.(( M.prinvars ./ sum(M.prinvars) ) * 100))
    return var_explained
end
function LifeHistoryEcoPCA(LH_df,EcoInfo,var="embryo")
    filteredEcoInfo = filter(x-> x.ID in LH_df.species, EcoInfo)
    basicecodf, fooddf, habitatdf = SimpleEcozoneInfo(filteredEcoInfo)
    # basicecodf = filter(x-> x.ID in lhidf.species, basicecodf)
    X = select(LH_df, Not(:species))
    # Perform PCA
    PCA = @load PCA pkg=MultivariateStats
    pca_model = PCA(maxoutdim=2)
    mach = machine(pca_model, X)
    fit!(mach)
    # Transform the data using the PCA model
    X_pca = MLJ.transform(mach, X)
    explained_vars = explained_var(Matrix(X))
    # Plot the first two components
    if var == "food"
        pcafig = scatter(X_pca[:, 1], X_pca[:, 2], group=fooddf[:,"adult"] ,xlabel="PC1 ($(explained_vars[1])%)", ylabel="PC2 ($(explained_vars[2])%)", title="Life History PCA, colored by $var",markerstrokewidth=0)

    elseif var == "habitat"
        pcafig = scatter(X_pca[:, 1], X_pca[:, 2], group=habitatdf[:,"adult"] ,xlabel="PC1 ($(explained_vars[1])%)", ylabel="PC2 ($(explained_vars[2])%)", title="Life History PCA, colored by $var",markerstrokewidth=0)
    else
        pcafig = scatter(X_pca[:, 1], X_pca[:, 2], group=basicecodf[:,var] ,xlabel="PC1 ($(explained_vars[1])%)", ylabel="PC2 ($(explained_vars[2])%)", title="Life History PCA, colored by $var",markerstrokewidth=0)
    end
    # pcafig = scatter(X_pca[:, 1], X_pca[:, 2], group=basicecodf[:,var] ,xlabel="PC1 ($(explained_vars[1])%)", ylabel="PC2 ($(explained_vars[2])%)", title="Life History PCA, colored by $var",markerstrokewidth=0)
    return pcafig
end
function LifeHistoryPCA2(LH_df,var="Wwi")
    X = Matrix(select(LH_df, Not(:species)))
    # Perform PCA
    M = fit(MultivariateStats.PCA,X, maxoutdim=2)
    # Transform the data using the PCA model
    # X_pca = MultivariateStats.predict(M, X)
    X_pca = M.proj
    # X_pca= MultivariateStats.reconstruct(M, X_pca)
    var_explained = Int.(round.(( M.prinvars ./ sum(M.prinvars) ) * 100))
    # Plot the first two components
    pcafig = scatter(X_pca[:, 1], X_pca[:, 2], marker_z=LH_df[:,var],color=:viridis ,xlabel="PC1 ($(var_explained[1])%)", ylabel="PC2 ($(var_explained[2])%)", title="Life History PCA, colored by $var", colorbar=true,markerstrokewidth=0)
    return pcafig
end
function LifeHistoryEcoPCA_customcolors(LH_df,EcoInfo,var="embryo")
    basicecodf, fooddf, habitatdf = SimpleEcozoneInfo(EcoInfo)
    basicecodf = filter(x-> x.ID in lhidf.species, basicecodf)
    X = select(LH_df, Not(:species))
    # Perform PCA
    PCA = @load PCA pkg=MultivariateStats
    pca_model = PCA(maxoutdim=2)
    mach = machine(pca_model, X)
    fit!(mach)
    # Transform the data using the PCA model
    X_pca = MLJ.transform(mach, X)
    explained_vars = explained_var(Matrix(X))
    # Plot the first two components
    colors = [:grey, :blue, :green]  # Define three distinct colors
    pcafig = scatter(
        X_pca[:, 1], X_pca[:, 2], 
        group=basicecodf[:, var], 
        xlabel="PC1 ($(explained_vars[1])%)", 
        ylabel="PC2 ($(explained_vars[2])%)", 
        title="Life History PCA, colored by $var",
        markerstrokewidth=0, 
        palette=colors  # Assign the custom color palette
    )
    return pcafig
end
#### Get the life history info ####
@time LH_df = GetLifeHistoryDFOld(AMPCOL)
@time LH_df_noguess = GetLifeHistoryDFOld(AMPCOL,true)
LH_df = LH_df_noguess
#### PCA ####
species_list = CSV.read(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\webscraping\recent_amp_scrape\AmP_species_list.csv", DataFrame)
filter!(x-> x.ID in unique(LH_df.species), species_list)
LifeHistoryPCA(LH_df,"Li")
LifeHistoryPCA(LH_df,"Ri")
ff, ddf = LifeHistoryPCA(LH_df,"Wwi")
LifeHistoryPCA(LH_df,"Wwb")
LifeHistoryPCA(LH_df,"ab")
ff1, ddf1 =LifeHistoryPCA(LH_df,"am")

# CSV.write(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\lucile_data\lifehistoryp_df.csv",LH_df )

# CSV.write(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\lucile_data\lifehistorypca.csv", DataFrame(ID=LH_df.species, PC1=ddf.x1, PC2=ddf.x2))

savefig("lifehistorypca.png")

#consider only Chordata
chordinds = species_list.Phylum .== "Chordata"
LH_df_chord = filter(x-> x.species in species_list.ID[chordinds], LH_df)
LifeHistoryPCA(LH_df_chord,"Li")
LifeHistoryPCA(LH_df_chord,"Ri")
LifeHistoryPCA(LH_df_chord,"Wwi")
LifeHistoryPCA(LH_df_chord,"Wwb")
LifeHistoryPCA(LH_df_chord,"ab")
LifeHistoryPCA(LH_df_chord,"am")

# This PCA Approach was done in doi:10.1038/s41559-019-0938-7
# Percentage of variance explained by each component should be added!!!
# NOTE: In order to do regression: Do standardisation and PCA on training data, then apply the same transformation to the test data
# life span, reproduction rate, weight at birth, ultimate weight, age at birth, and ultimate length

########### Regression ###########
mito_df = CSV.read(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\webscraping\recent_amp_scrape\AMPMITO.csv", DataFrame)
yt,yp = MitoGenomeLHIRegression(mito_df,LH_df,kfolds=50,kmers=4,model=rego)
r2 = rsquared(Float64.(yt), Float64.(yp))
println("R-squared: $r2")

# new more customized version
@time lhidf2 = GetLifeHistoryDF(AMPCOL,vars=["am","Wwi","Ri","Li"],removeguess=true) #15.528210 seconds
LifeHistoryPCA(lhidf2,"Li")

### use MVS to do the PCA (and include % of variance explained) ###
using MultivariateStats
X = select(lhidf2, Not(:species))
pca_model = PCA(maxoutdim=2)
mach = machine(pca_model, X)

M = fit(MultivariateStats.PCA, Matrix(X))
M.prinvars
var_explained = Int.(round.(( M.prinvars ./ sum(M.prinvars) ) * 100))
M.tvar
M.proj

scatter(M.proj[:, 1], M.proj[:, 2], marker_z=lhidf2[!,:Li],color=:viridis ,xlabel="PC1 ($(var_explained[1])%)", ylabel="PC2 ($(var_explained[2])%)", title="Life History PCA, colored by Li", colorbar=true,markerstrokewidth=0)


 # Get the loadings (eigenvectors)
 loadings = M.prinvars
 # Plot the first two components
 pcafig = scatter(X_pca[:, 1], X_pca[:, 2], marker_z=LH_df[:, var], color=:viridis, xlabel="PC1", ylabel="PC2", title="Life History PCA, colored by $var", colorbar=true, markerstrokewidth=0)
 # Add directional vectors
 for i in 1:size(loadings, 1)
     arrow!(pcafig, [0, 0], [loadings[i, 1], loadings[i, 2]], label=names(X)[i], color=:black, linewidth=2, arrow=:arrow)
 end

### consider ecozone info ###
EcoInfo = CSV.read(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\webscraping\recent_amp_scrape\AmP_EcoInfo.csv", DataFrame)
basicecodf, fooddf, habitatdf = SimpleEcozoneInfo(EcoInfo)

EcoInfo_limited = filter(x-> x.ID in LH_df.species, EcoInfo)
basicecodf_limited = filter(x-> x.ID in LH_df.species, basicecodf)
fooddf_limited = filter(x-> x.ID in LH_df.species, fooddf)
habitatdf_limited = filter(x-> x.ID in LH_df.species, habitatdf)
species_list_limited = filter(x-> x.ID in LH_df.species, species_list)

CSV.write(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\lucile_data\ecoinfo.csv", EcoInfo_limited)
CSV.write(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\lucile_data\basicecodf.csv", basicecodf_limited)
CSV.write(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\lucile_data\fooddf.csv", fooddf_limited)
CSV.write(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\lucile_data\species.csv", species_list_limited)
CSV.write(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\lucile_data\habitatdf.csv", habitatdf_limited)


##### Ecozone annotated PCA ####
display_table(species_list)
lhidf = GetLifeHistoryDF(AMPCOL) 

basicecodf2 = filter(x-> x.ID in lhidf.species, basicecodf)

M = fit(MultivariateStats.PCA, Matrix(select(lhidf, Not(:species))))

M.prinvars
M.tvar

PC1 = M.proj[:,1]
PC2 = M.proj[:,2]
PC3 = M.proj[:,3]

scatter(PC1, PC2,PC3, group=basicecodf2.embryo, legend=:topleft)
scatter(PC1, PC2, group=basicecodf2.gender, legend=:topleft)
scatter(PC1, PC2, group=basicecodf2.reproduction, legend=:topleft)
scatter(PC1, PC2, group=basicecodf2.ecozone, legend=:topleft)
scatter(PC1, PC2, marker_z=lhidf[:,"Li"],color=:viridis)


LifeHistoryPCA2(lhidf,"Ri")
LifeHistoryPCA2(lhidf,"ab")
LifeHistoryPCA2(lhidf,"am")
pca_model = PCA(maxoutdim=2)
mach = machine(pca_model, X)

M.prinvars
LifeHistoryEcoPCA(LH_df,EcoInfo,"habitat")






figo,xpca = LifeHistoryPCA(lhidf)
p1 = scatter(xpca.x1, zeros(length(xpca.x1)),size=(800, 100) , marker_z=lhidf[:,"am"],c=:viridis, markersize=5, label="am", xlabel="PCA1", ylabel="", yticks=:none, ylims=(-1,1), markerstrokewidth=0)
p2 = scatter(xpca.x1, zeros(length(xpca.x1)),size=(800, 100) , marker_z=lhidf[:,"ab"],c=:plasma, markersize=5, label="ab", xlabel="PCA1", ylabel="", yticks=:none, ylims=(-1,1), markerstrokewidth=0)
p3 = scatter(xpca.x1, zeros(length(xpca.x1)),size=(800, 100) , marker_z=lhidf[:,"Ri"],c=:coolwarm, markersize=5, label="Ri", xlabel="PCA1", ylabel="", yticks=:none, ylims=(-1,1), markerstrokewidth=0)
p4 = scatter(xpca.x1, zeros(length(xpca.x1)),size=(800, 100) , marker_z=lhidf[:,"Wwi"],c=:PuOr, markersize=5, label="Wwi", xlabel="PCA1", ylabel="", yticks=:none, ylims=(-1,1), markerstrokewidth=0)
p5 = scatter(xpca.x1, zeros(length(xpca.x1)),size=(800, 100) , marker_z=lhidf[:,"Wwb"], markersize=5, label="Wwb", xlabel="PCA1", ylabel="", yticks=:none, ylims=(-1,1), markerstrokewidth=0)
p6 = scatter(xpca.x1, zeros(length(xpca.x1)),size=(800, 100) , marker_z=lhidf[:,"Li"], markersize=5, label="Li", xlabel="PCA1", ylabel="", yticks=:none, ylims=(-1,1), markerstrokewidth=0)

plot(p1, p2, p3, p4, p5, p6, layout=@layout([a; b; c; d; e; f]), size=(800, 600))



bb,ff,hh = SimpleEcozoneInfo(EcoInfo)
display_table(GetFoodTable(EcoInfo))
