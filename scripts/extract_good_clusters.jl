push!(LOAD_PATH, dirname(Base.PROGRAM_FILE)*"/julia_modules/")
using sanitise
using ArgParse2
using DataFrames
using CSV 
function julia_main()::Cint
    parser = ArgumentParser( description = "Extract clusters of a proper size and without close chimeras")
    #add_argument!(parser, "-i", "--input-fasta", type = String, help = "Input file")
    add_argument!(parser, "-a", "--cluster-lower-level", type = String, help = "Clustering to group reads for assembly")
    add_argument!(parser, "-b", "--cluster-higher-level", type = String, help = "Clustering to group reads to check 4 chimeras")
    add_argument!(parser, "-t", "--tax-assign", type = String, help = "Taxonomy assignment files - down to genus level. Should match the lower level clustering fasta file ")
    add_argument!(parser, "-m", "--min-cluster-size", type = Float64, help = "Minimum cluster size")
    add_argument!(parser, "-f", "--fasta", type = String, help = "Fasta file matching lower lever of clustering ")
    add_argument!(parser, "-c", "--detected-chimers", type = String, help = "Info of detected chimeric clusters ", default = "chimeric_info.csv")
    args = parse_args(parser)
    aln_all = Alignment(args.fasta) 
    centroid_genus = parse_dada_results(args.tax_assign)
    cl1  = read_jc_clusters(args.cluster_lower_level)
    #check which cluster has proper size
    proper_sizes = mark_proper_size_clusters(cl1,args.min_cluster_size)
    #get higher level cluster - it should group such close species as E. Coli & Salmonela and their chimeric reads
    cl2   = read_jc_clusters(args.cluster_higher_level)
    recomb_data = DataFrame()
    groupct = 0
    for cl in keys(cl2) 
        ids_for_common_analyses = Array{String,1}()
        for rep in cl2[cl]
                if proper_sizes[rep]
                    push!(ids_for_common_analyses,rep)
                end 
        end 
        if length(ids_for_common_analyses) > 2
            #a case where chimeric cluster can be formed
            groupct += 1
            #println("------------------")
            #for id in ids_for_common_analyses 
            #    println(id," ",length(cl1[id])," ",centroid_genus[id]," ", proper_sizes[id])
            #end 
            sub_aln = sub_alignment(aln_all,ids_for_common_analyses)
            recombs = eval_recomb(sub_aln)
            recombs[:Group] = fill(groupct,nrow(recombs))
            recomb_data = vcat(recomb_data,recombs)
        end 
    end
    println(stderr,"Chimeric reads detection info: ",args.detected_chimers)
    CSV.write(args.detected_chimers, recomb_data)

    
    #aln = sanitise.Alignment(args.input_fasta)
    #println(aln_conservation(aln,ignore_amounts=true))

    return 0
end

julia_main()
#mycoolfunction()
