function go()
    f2 = ARGS[1] #names of sorted by size cenroids of clusters
    println(stderr,"Reading cluster mapping from $f2")
    map =Dict{String,Array{String,1}}()
    cur=""
    for l in readlines(f2)
        println(l)
        parts = split(l)
        if parts[2] !=cur
            map[parts[2]]=Array{String,1}()
            cur = parts[2]
        end
        push!(map[parts[2]],parts[1])
    end
    sizes = Array{Tuple{String,Int64},1}()
    for k in keys(map)
        push!(sizes,(k,length(map[k])))
    end
    sort!(sizes,by=x->x[2],rev=true)
    for p in map[sizes[1][1]]
        println(p)
    end
end
go()
