using ArgParse2
using DataFrames 
using DataFramesMeta 
using CSV 
using Statistics 

function main()
    parser = ArgumentParser(prog = "",
                        description = "Evaluates ckstering accuracy based on expected gebus distribution",
                        epilog = "sizes..ajusted")
    #add_argument!(parser, "--read-counts","-c", help = "directory with Reads counts based on aligned clustering",type = String)
    add_argument!(parser, "--reference-values","-r", help = "Reference values",type = String)
    add_argument!(parser,"--read-counts","-c", metavar="COUNTS", nargs = "+", required = true)
    args = parse_args(parser)
    ref_data = CSV.read(args.reference_values)
    ct = 0 
    for cntf in args.read_counts        
        ct += 1
        print(stderr,ct," ",cntf,"\n") 
        counts = DataFrame(CSV.read(cntf)) 
        counts = @transform( counts, Genus2 = getindex.(split.(:Genus,"/"),1) )
        counts = @select(counts, :Genus2, :Maximum_abundance)
        dd =join(counts,ref_data, on = :Genus2 => :Genus, kind = :inner)
        coro = cor( dd[:,:Maximum_abundance]  ,  dd[:,Symbol("16S")] )
        println(cntf," ",coro," ",nrow(dd))
    end 
end

main()

