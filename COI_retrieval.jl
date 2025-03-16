using CSV, DataFrames

function reorganise_coi_df(coi_DF)
    # coi_DF[!,"identifier"] = coi_DF.accession .* " " .* coi_DF.species
    coi_df = select(coi_DF, [:identifier, :sequence])
    specieses = []
    accessions = ["" for i in 1:size(coi_df)[1]]
    samples = []
    rest = []
    genes = []
    for i in 1:size(coi_df)[1]
        splits_ =  split(coi_df.identifier[i], '|')
        sample,species,gene = splits_[1:3]
        push!(specieses,strip(species))
        push!(samples,strip(sample))
        push!(genes,strip(gene))
        if length(splits_) == 4
            accessions[i] = strip(splits_[4])
        end
    end
    coi_df[!,"species"] = specieses
    coi_df[!,"accession"] = accessions
    coi_df[!,"sample"] = samples
    coi_df[!,"gene"] = genes
    return select(coi_df, [:sample,:species,:gene,:accession, :sequence])
end

function get_coi(path)
    coi_df = parse_fasta_to_dataframe(path)
    coi_df[!,"identifier"] = coi_df.accession .* " " .* coi_df.species
    coi_df = select(coi_df, [:identifier, :sequence])
    return reorganise_coi_df(coi_df)
end
function sequence_completeness(sequence)
    Ncount = count(c -> c == 'N', sequence)
    dashcount = count(c -> c == '-', sequence)
    return 1 - (Ncount + dashcount)/length(sequence)
end
function remove_identical_seqs(fasta_df::DataFrame)
    # Find the unique sequences and their first occurrences
    unique_sequences = unique(fasta_df.sequence)
    first_occurrences = Dict(seq => findfirst(x -> x == seq, fasta_df.sequence) for seq in unique_sequences)
    
    # Create a new DataFrame with only the first occurrences of each unique sequence
    unique_df = fasta_df[first.(values(first_occurrences)), :]
    
    return unique_df
end
### Loading datasets ###

scraped_amp = CSV.read(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\webscraping\recent_amp_scrape\AmP_species_list.csv", DataFrame)
@time chordata_COI_df = parse_fasta_to_dataframe(raw"C:\Users\msvhaeve\Downloads\chordataCOI.fas")
chordata_COI_df[!,"identifier"] = chordata_COI_df.accession .* " " .* chordata_COI_df.species
chordata_COI_df = select(chordata_COI_df, [:identifier, :sequence])
malacopath = raw"C:\Users\msvhaeve\Downloads\Malacostraca.fas"
branchiopath = raw"C:\Users\msvhaeve\Downloads\BranchiopodaCOI.fas"
echinopath = raw"C:\Users\msvhaeve\Downloads\EchinodermataCOI.fas"
molluscopath = raw"C:\Users\msvhaeve\Downloads\MolluscaCOI.fas"
cnidarpath = raw"C:\Users\msvhaeve\Downloads\cnidarioCIO.fas"
COI_table_amp = CSV.read(raw"C:\Users\msvhaeve\Documents\DEBpaper\DEBPaper_datasets\COI_df_l10_sl1000_annotated.csv", DataFrame) #previous database

### testing grounds ###


@time Chordata_COI_df = reorganise_coi_df(chordata_COI_df)
Chordata_COI = filter(x -> x.gene == "COI-5P", Chordata_COI_df) 
@time Chordata_COI_AMP = get_amp_mitototal_df(Chordata_COI,scraped_amp) #242 seconds

genecount = combine(groupby(Chordata_COI_df, :gene), nrow => :Count)
display_table(genecount)


Chordata_COI_AMP[!,"Completeness"] = sequence_completeness.(Chordata_COI_AMP.sequence)
Chordata_COI_AMP[!,"Length"] = length.(Chordata_COI_AMP.sequence)

completenesses = sort(Chordata_COI_AMP.Completeness, rev=true)
lengths = sort(Chordata_COI_AMP.Length, rev=true)

# using Plots
# histogram(completenesses, bins=100, xlabel="Sequence completeness", ylabel="Number of sequences", title="Sequence completeness of Chordata COI sequences")
# histogram(lengths, bins=100, xlabel="Sequence length", ylabel="Number of sequences", title="Sequence length of Chordata COI sequences", yaxis=:log)

Chordata_COI_AMP2 = filter(x -> x.Completeness > 0.60, Chordata_COI_AMP)
Chordata_COI_AMP3 = filter(x -> 500 < x.Length < 800, Chordata_COI_AMP2)
Chordata_COI_AMP4 = filter(x -> 600 < x.Length < 750, Chordata_COI_AMP2)

@time Chordata_COI_AMP0 = remove_identical_seqs(Chordata_COI_AMP)
unique(Chordata_COI_AMP3.ID)
Chordata_COI_AMP2 = filter(x -> x.Completeness > 0.60, Chordata_COI_AMP0)
display_table(Chordata_COI_AMP)
maximum(Chordata_COI_AMP.Length)

malaco = get_coi(malacopath)
branchio = get_coi(branchiopath)
echino = get_coi(echinopath)
mollusco = get_coi(molluscopath)
cnidario = get_coi(cnidarpath)

nonchords = vcat(malaco,branchio,echino,mollusco,cnidario)
nonchordsCOI = filter(x -> x.gene == "COI-5P", nonchords)

@time nonchordsCOI_AMP = get_amp_mitototal_df(nonchordsCOI,scraped_amp) #242 seconds
@time nonchordsCOI_AMP0 = remove_identical_seqs(nonchordsCOI_AMP)
nonchordsCOI_AMP0[!,"Completeness"] = sequence_completeness.(nonchordsCOI_AMP0.sequence)
nonchordsCOI_AMP0[!,"Length"] = length.(nonchordsCOI_AMP0.sequence)
AMP_COI = vcat(Chordata_COI_AMP0,nonchordsCOI_AMP0)

# CSV.write(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\webscraping\recent_amp_scrape\Bold_AMP_COI.csv", AMP_COI)

# Combine with previous COI data (from COInr database)
prev_coi = get_amp_mitototal_df(select(COI_table_amp,[:species,:seqID,:length,:sequence]),scraped_amp)
extra_coi = filter(x->x.ID in setdiff(unique(prev_coi.ID),unique(AMP_COI.ID)),prev_coi)
extra_coi0 = remove_identical_seqs(extra_coi)
extra_coi0[!,"Completeness"] = sequence_completeness.(extra_coi0.sequence)
extra_coi0[!,"Length"] = length.(extra_coi0.sequence)
setdiff(names(AMP_COI), names(extra_coi0))

extra_coi0[!,"accession"] = extra_coi0.seqID
extra_coi0[!,"gene"] .= "COI-5P"
extra_coi0[!,"sample"] .= ""

Bold_AMP_COI_Final = vcat(AMP_COI,select(extra_coi0,Symbol.(names(AMP_COI))))

# CSV.write(raw"C:\Users\msvhaeve\OneDrive - UGent\Desktop\webscraping\recent_amp_scrape\Bold_AMP_COI_Final.csv", Bold_AMP_COI_Final)