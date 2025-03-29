##### Obtain a mitogenome dataset using different sources #####

### Load packages and define relevant functions ###
using CSV, DataFrames

display_table(df) = VSCodeServer.vscodedisplay(df)
function parse_fasta_to_dataframe(fasta_file::String)
    # Initialize an empty DataFrame
    df = DataFrame(accession = String[], species = String[], sequence = String[])
    
    # Initialize variables to store the current record
    accession = ""
    species = ""
    sequence = ""
    
    # Open and read the FASTA file
    open(fasta_file, "r") do file
        for line in eachline(file)
            if startswith(line, ">")
                # If we encounter a new header, save the previous record
                if accession != ""
                    push!(df, (accession, species, sequence))
                end
                
                # Extract accession and species from the header
                # Assuming header format: ">accession species"
                parts = split(line[2:end], ' ')
                accession = parts[1]
                species = join(parts[2:end], ' ')
                
                # Reset sequence for the new record
                sequence = ""
            else
                # Append the line to the sequence
                sequence *= strip(line)
            end
        end
        
        # Don't forget to add the last record
        if accession != ""
            push!(df, (accession, species, sequence))
        end
    end
    
    return df
end
function get_mitofish_df(mitofish_dir)
    mitofish_files = sort(readdir(mitofish_dir))
    fish_files = joinpath.(mitofish_dir, mitofish_files)
    df_mitos = parse_fasta_to_dataframe(fish_files[1])
    a,gi_,r, ref_, t, tax_, name_ = split(df_mitos[1,1],'|')
    gis = [gi_]
    refs = [ref_]
    species = [name_*" "*df_mitos[1,2]]
    for file in fish_files[2:end]
        df = parse_fasta_to_dataframe(file)
        a,gi,r,ref,t,tax,name = split(df[1,1],'|')
        push!(gis, gi)
        push!(refs, ref)
        push!(species, name*" "*df[1,2])
        df_mitos = vcat(df_mitos, parse_fasta_to_dataframe(file))
    end
    df_mitos[!, :gi] = gis
    df_mitos[!, :ref] = refs
    df_mitos[!, :species] = species
    df_mitos[!,:accession] = df_mitos.ref
    df_mitos[!,:ID] = [split(x, ' ')[1]*'_'*split(x, ' ')[2] for x in df_mitos.species] 
    return select(df_mitos, [:accession, :species, :sequence])
end
function get_name(file)
    namesplit = split(file[1:end-3],'_')[3:end]
    return namesplit[1]*'_'*namesplit[2]
end
# Function to process species name
function process_species_name(name)
    name = replace(name, r"mitochondrion" => "")
    name = replace(name, r"complete" => "")
    name = replace(name, r"genome" => "")
    name = replace(name, r"voucher.*" => "")
    name = replace(name, r"sequence" => "")
    name = replace(name, r"mitogenome" => "")
    name = replace(name, r"mitochondrial" => "")
    name = replace(name, r"genomic" => "")
    name = replace(name, r"nearly" => "")
    name = replace(name, r"isolate.*" => "")
    name = replace(name, r"haplotype.*" => "")
    name = replace(name, r"DNA*" => "")
    name = replace(name, r"assembly.*" => "")
    name = replace(name, r"organelle.*" => "")
    name = replace(name, r" x .*" => "")
    name = replace(name, r" [A-Z] " => " ")
    name = replace(name, r"strain.*" => "")
    name = replace(name, r"\d+" => "")  # Remove all numbers
    name = replace(name, r"[[:punct:]]" => "")  # Remove all punctuation

    return strip(name)
end

is_NC(value) = occursin(r"NC_", value)
function process_crosses(name)
    name = replace(name, r" x .*" => " ")
    return strip(name)
end
function accession_enrichment(AMPMITO)
    NC_AMPMITO = filter(x -> is_NC(x.accession), AMPMITO)
    XC_AMPMITO = filter(x -> !is_NC(x.accession), AMPMITO)
    unique_non_NC = filter(x -> x.ID in setdiff(XC_AMPMITO.ID,NC_AMPMITO.ID), XC_AMPMITO)
    return vcat(NC_AMPMITO,unique_non_NC)
end
function ID_unification(AMPMITO)
    unique_IDs = unique(AMPMITO.ID)
    retained_inds = []
    for ID in unique_IDs
        push!(retained_inds,findfirst(occursin.(AMPMITO.ID, ID)))
    end
    return AMPMITO[retained_inds,:]
end
function correct_accession(name)
    name = replace(name, r"(\.1).*" => s"\1")  # Remove everything after ".1"
    name = replace(name, r"(\.2).*" => s"\1") 
    return name
end
function get_amp_mitototal_df(mitototal,scraped_amp)
    mitototal.species = process_crosses.(mitototal.species)
    IDs_amp = []
    corresponding_entries = []
    for i in 1:size(mitototal)[1]
        push!(IDs_amp,findfirst(occursin.(scraped_amp.ScientificName, mitototal.species[i])))
        push!(corresponding_entries,i)
    end
    AMP_inds = IDs_amp[.!isnothing.(IDs_amp)]
    Mito_inds = corresponding_entries[.!isnothing.(IDs_amp)]
    amp_mitototal = scraped_amp[AMP_inds,:]
    mitos_mitototal = mitototal[Mito_inds,:]
    AMPMITO = hcat(amp_mitototal,mitos_mitototal)
    AMPMITO.species = process_species_name.(AMPMITO.species)
    return AMPMITO
end
function extend_amp_mitototal_df(amp_mitototal,new_mito_df,scraped_amp)
    new_mito_df = select(new_mito_df, [:accession, :species, :sequence])
    new_mito_df.species = process_crosses.(new_mito_df.species)
    new_mitos_df = get_amp_mitototal_df(new_mito_df,scraped_amp)
    extramito = filter(x -> x.accession in setdiff(new_mitos_df.accession,amp_mitototal.accession) , new_mitos_df)
    extramito.species = process_species_name.(extramito.species)
    extramito.accession = correct_accession.(extramito.accession)
    AMPMITO = vcat(amp_mitototal,extramito)
    return AMPMITO
end

#### Get Data #### 
## Load AmP data and mitogenomes downloaded from NCBI, MitoFish, and MitoAge 
include("SystemSpecific_filepaths.jl")
#### Construct total mitogenome dataset from the different sources #### 83.281899 seconds
@time begin
AMPMITO = get_amp_mitototal_df(MitoFish,scraped_amp)
AMPMITO = extend_amp_mitototal_df(AMPMITO,fastaa,scraped_amp)
AMPMITO = extend_amp_mitototal_df(AMPMITO,mitototal,scraped_amp)
AMPMITO = extend_amp_mitototal_df(AMPMITO,mammals,scraped_amp)
AMPMITO = extend_amp_mitototal_df(AMPMITO,birds,scraped_amp)
AMPMITO = extend_amp_mitototal_df(AMPMITO,prev_mitos,scraped_amp)
AMPMITO = extend_amp_mitototal_df(AMPMITO,allmito,scraped_amp)
AMPMITO = extend_amp_mitototal_df(AMPMITO,totalmito,scraped_amp)
AMPMITO = extend_amp_mitototal_df(AMPMITO,formermito,scraped_amp)
AMPMITO = extend_amp_mitototal_df(AMPMITO,MitoAge_df,scraped_amp)
AMPMITO = accession_enrichment(AMPMITO)
AMPMITO = ID_unification(AMPMITO)
end
