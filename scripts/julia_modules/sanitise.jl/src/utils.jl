"Readd clustering from jc/gc  file"
function read_jc_clusters(inf::String)
    map =OrderedDict{String,Array{String,1}}()
    pairs = Array{Tuple{String,String},1}()
    for l in readlines(inf)
	parts = split(l)
	push!(pairs,(parts[1],parts[2]))
    end
    #sort to prevent constant check if key exists
    sort!(pairs, by = x -> x[2])
    cur=""
    for parts in pairs 
	if parts[2] !=cur
	    map[parts[2]]=Array{String,1}()
	    cur = parts[2]
	end
	push!(map[parts[2]],parts[1])
    end
    #sort by the size of cluster
    sort!(map, by= x -> length(map[x]),rev=true)

    return(map)
end



"Return a bool array which length is equal to input OrderedDict keys and false mark small cluster"
function mark_proper_size_clusters(cls::OrderedDict{String,Array{String,1}},minfrac::Float64)
    cls_sizes = map(length,values(cls))
    all_ct = sum(cls_sizes)
    fracs = (cls_sizes ./ all_ct)  .>= minfrac 
    out = OrderedDict(zip(keys(cls),fracs))
    return(out)
end 

"Parase genus assignment from a dada results file"
function parse_dada_results(fd::String)
    out = 
    centroid_genus_map = Dict{String,String}()
    ct = 0
    nact = 0
    totalnr = 0
    genuses = Set{String}()
    genuses_counts = Dict{String,Int64}()
    for l in readlines(fd)
	ct += 1
	if ct > 1
	    parts = split(replace(l,"\""=>""))
	    rn = split(parts[1],";")[1]
	    size = parse(Int64,split(split(parts[1],";")[2],"=")[2])
	    totalnr += size
	    gn = parts[7]
	    if gn == "NA"
		nact += 1
		gn = "$gn"*"_"*"$nact"
	    else 
		if gn in genuses 
		    genuses_counts[gn] += 1
		else 
		    push!(genuses,gn)
		    genuses_counts[gn] = 1
		end 
		gnnr = genuses_counts[gn]
		gn = "$gn"*"_"*"$gnnr"
	    end
	    centroid_genus_map[rn] = gn
	end
    end
    return(centroid_genus_map)
end 



#
