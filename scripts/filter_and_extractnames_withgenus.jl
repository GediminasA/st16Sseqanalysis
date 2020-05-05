using DataStructures
function go()
f1 = ARGS[1] #names of sorted by size cenroids of clusters
f2 = ARGS[2] #cluster mapping file sorted by the secod column
f3 = ARGS[3] #dada2 assignment
outdir = ARGS[4] #directory for writing out clusters
minimumclustersize = parse(Float64,ARGS[5]) #directory for writing out clusters
println(stderr,"Reading cluster mapping from $f3")
centroid_genus_map = OrderedDict{String,Array{Any,1}}()
ct = 0
nact = 0
totalnr = 0
for l in readlines(f3)
    ct += 1
    if ct > 1
        parts = split(replace(l,"\""=>""))
        rn = split(parts[1],";")[1]
        size = parse(Int64,split(split(parts[1],";")[2],"=")[2])
        totalnr += size
        gn = parts[7]
        if gn == "NA"
            nact += 1
            gn = "$gn$nact"
        end
        if !(gn in keys(centroid_genus_map))
            centroid_genus_map[gn] = [Array{String,1}(),0]
        end
        push!(centroid_genus_map[gn][1],rn)
        centroid_genus_map[gn][2] += size
    end
end
#sort the dictionary by size
centroid_genus_map_array = Array{Array{Any,1},1}()
for v in keys(centroid_genus_map)
    push!(centroid_genus_map_array,[v,centroid_genus_map[v][1],centroid_genus_map[v][2]])
end

sort!(centroid_genus_map_array,by=x->x[3],rev=true)
#for i in centroid_genus_map_array
#    println(i[1]," ",i[3])
#end



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
finfo = open("$outdir/cluster_genus_size.csv","w")
println(finfo,"Cluster;Genus;Size")
for genusinf in centroid_genus_map_array
    ct += 1
    gn = replace(genusinf[1],"/"=>"-")
    size = genusinf[3]
    if size >= minimumclustersize*totalnr
        out = "$outdir/$ct"
        f=open(out,"w")
        for centr in genusinf[2]
            for id in map[centr]
                println(f,id)
            end
        end
        close(f)
        f=open("$outdir/$ct.$gn","w")
        println(f,"$size")
        close(f)
        println(finfo,"$ct;$gn;$size")
    else
        println("Skiping cluster $ct $gn as too small size less than $minimumclustersize fraction")
    end
end
close(finfo)
end
go()
