push!(LOAD_PATH, "scripts/julia_modules/ArgParse2.jl/src")
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
    ref_data =  DataFrame!(CSV.File(args.reference_values))
    println("File CorvsMaximumcluster CorvsTotal CorvsNumberOfClusters NumberOfClusters DetectedReqGenus")
    ct = 0
    cor1 = 0 
    cor2 = 0 
    cor3 = 0 
    for cntf in args.read_counts        
        ct += 1
        counts = DataFrame!(CSV.File(cntf)) 
        counts = @transform( counts, Genus2 = getindex.(split.(:Genus,"/"),1) )
        counts = @transform( counts, Genus2 = getindex.(split.(:Genus2,"-"),1) )
        counts = @transform( counts, Genus2 = getindex.(split.(:Genus2),1) )
        counts1 = @select(counts, :Genus2, :Maximum_abundance)
        counts2 = @select(counts, :Genus2, :Total_abundance)
        counts3 = @select(counts, :Genus2, :Number_of_sequences)
        dd1 =join(counts1,ref_data, on = :Genus2 => :Genus, kind = :inner)
        dd2 =join(counts2,ref_data, on = :Genus2 => :Genus, kind = :inner)
        dd3 =join(counts3,ref_data, on = :Genus2 => :Genus, kind = :inner)
        if (nrow(dd1) > 1)
            cor1 = cor( dd1[:,:Maximum_abundance]  ,  dd1[:,Symbol("16S")] )
            cor2 = cor( dd2[:,:Total_abundance]  ,  dd2[:,Symbol("16S")] )
            cor3 = cor( dd3[:,:Number_of_sequences]  ,  dd3[:,Symbol("16S")] )
        end 
        println(cntf," ",cor1," ",cor2," ",cor3," ", sum(dd3[:,:Number_of_sequences])," ",nrow(dd1))
    end 
end

main()

