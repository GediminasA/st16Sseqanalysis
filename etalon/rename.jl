function run()
    out="Zymo1-2X-65C_S6_contigs.fasta"
    ref="rRNA_ampl_1500Genomic.fasta"
    mapf = "map"
    #read names mapping
    text=read(open(ref,"r"),String)
    f=open(mapf,"r")
    ct = 0
    for l in readlines(mapf)
        o,r = split(l)
        text = replace(text,o=>r)

    end

    outf=open(out,"w")
    close(outf)
end

run()
