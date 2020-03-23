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
ct = 0
for l in readlines(f1)
    ct +=1
    org_name = split(l,";")[1]
    cluster = map_new[org_name]
    println("$org_name\t$cluster")
end
println(stderr, "$ct cluster were written out")
end
go()
