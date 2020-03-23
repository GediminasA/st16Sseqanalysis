function go()
f1 = ARGS[1]
f2 = ARGS[2]
println(stderr,"Reading new cluster mapping from $f2")
map_new =Dict{String,String}()
for l in readlines(f2)
    parts = split(l,"\t")
    map_new[split(parts[1],";")[1]]=split(parts[2],";")[1]

end
println("Readed in nea mapping")
skipped = 0
for l in readlines(f1)
    parts = split(l,"\t")
    org_name = split(parts[1],";")[1]
    old_cluster = split(parts[2],";")[1]
    new_cluster = get(map_new, old_cluster, nothing)
    if new_cluster == nothing
        skipped += 1
    else
        println("$org_name\t$new_cluster")
    end
end
println(stderr, "$skipped cluster were discarded as were not found in the new mapping")
end
go()
