using DataFrames
using DataFramesMeta
using CSV

blastfile = ARGS[1]
outstem = ARGS[2]

blast_fields ="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
df = CSV.read(blastfile,header = Array{String,1}(split(blast_fields)))
sizes = fill(0,size(df)[1])
species = fill("",size(df)[1])
read_name = fill("",size(df)[1])
for i in 1:size(df)[1]
    sizes[i] = parse(Int64,split(df.qseqid[i],"=")[2])
    read_name[i] = split(df.qseqid[i],";")[1]
    species[i] = String(split(df.sseqid[i],":")[1])
end
df.sizes = sizes
df.species = species
df.read_name = read_name

dfout = DataFrame(Species = Array{String,1}(),
    Maximum_length = Array{Int64,1}(),
    Number_of_sequences = Array{Int64,1}(),
    Maximum_abundance = Array{Int64,1}(), #parses size argument in name, not the real number od sequences
    Total_abundance = Array{Int64,1}(), #parses size argument in name, not the real number od sequences
    Samples = Array{String,1}())
for s in unique(species)
    dfs = @linq df |>
               where(:species .== s) |>
               orderby(-:sizes)
    dfs_ids = @linq dfs |>
        select(:read_name)
    out = outstem*"_ids_of_"*s
    CSV.write(out,dfs_ids,header=[""])
    maxl = maximum(dfs.qlen)
    maxnb = maximum(dfs.sizes)
    totalnb = sum(dfs.sizes)
    smp = split(outstem,'/')[end]
    push!(dfout,[s,maxl,size(dfs)[1],maxnb, totalnb, smp])
end
out = outstem*"_blast_summary.tsv"
CSV.write(out,dfout,delim='\t')
