ucinfile  = ARGS[1]
curr =" "
for l in readlines(ucinfile)
    if l[1] == 'H'
        parts = split(l)
        mem,rep = split(parts[9],";")[1],split(parts[10],";")[1]
        println("$mem\t$rep")
    end
    if l[1] == 'C'
        parts = split(l)
        mem = split(parts[9],";")[1]
        println("$mem\t$mem")
    end
end
