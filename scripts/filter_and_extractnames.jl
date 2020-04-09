function go()
f1 = ARGS[1] #names of sorted by size cenroids of clusters
f2 = ARGS[2] #cluster mapping file sorted by the secod column
outdir = ARGS[3] #directory for writing out clusters

println(stderr,"Reading cluster mapping from $f2")
map =Dict{String,Array{String,1}}()
cur=""
for l in readlines(f2)
    parts = split(l)
    if parts[2] !=cur
        map[parts[2]]=Array{String,1}()
        cur = parts[2]
    end
    push!(map[parts[2]],parts[1])
end
nb_of_cl = length(keys(map))
println("Readed in maping for $nb_of_cl number of clusters")
println("File names for clusters will be writen out to the directory $outdir")
ct = 0
for l in readlines(f1)
    nm = split(l)[1]
    ct += 1
    out = "$outdir/$ct"
    f=open(out,"w")
    for id in map[nm]
        println(f,id)
    end
    close(f)
end
end
go()
