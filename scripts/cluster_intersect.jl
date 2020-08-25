using BioSequences
using CodecZlib
using FastaIO
using FASTX 
using Base.Threads
using ProgressMeter
push!(LOAD_PATH, "scripts/julia_modules/ArgParse2.jl/src")
println()
using ArgParse2



function write_out(centroids::Array{String,1},mapping::Array{Tuple{Array{String,1},String},1},outf::String)
    fo = open(outf,"w")
    println("Riting out to $outf  ...")
    @showprogress 1 "Writing out to $outf ..." for i in 1:length(centroids)
    #for i in 1:length(centroids)
        seq = centroids[i]
        name = mapping[i][2]
        print(fo,">$name\n$seq\n")
    end
    close(fo)
end

function write_out_fastq(seqs::Dict{String,FASTX.FASTQ.Record}, mapping::Array{Tuple{Array{String,1},String},1},outf::String)
    fo = FASTQ.Writer(GzipCompressorStream(open(outf,"w")))
    println("Riting out to $outf  ...")
    @showprogress 1 "Writing out to $outf ..." for pair in mapping 
   # for pair in mapping 
        centroid = seqs[pair[2]]
        for seq_id in pair[1]
            cluster_mem =  FASTQ.Record(seq_id, FASTX.FASTQ.sequence(centroid), FASTX.FASTQ.quality(centroid))
            write(fo,cluster_mem)
        end 
        #seq = centroids[i]
        #name = mapping[i][2]
        #print(fo,">$name\n$seq\n")
    end
    close(fo)
end

function run_clustering(seqs::Dict{String,FASTX.FASTQ.Record}, mapping::Array{Array{String,1},1})
    seqs2 = Dict{String,String}()
    for k in keys(seqs)
        seqs2[k] = FASTX.FASTQ.sequence(seqs[k])
    end 
    run_clustering(seqs2,mapping)
end
 
function run_clustering(seqs::Dict{String,String}, mapping::Array{Array{String,1},1})
    #mapping = mapping[1:10000]
    n = length(mapping)
    p = Progress(n)
    println(stderr,"Starting $n clustering runs...")
    centroids = Array{String,1}(undef, n)
    Threads.@threads for i in 1:n
                r1s = mapping[i]
                if length(r1s) == 1
                    centroids[i] = seqs[r1s[1]]
                else
                    centroids[i] = cluster_r2(r1s, seqs,t=3)
                end
                next!(p)
    end
    return(centroids)
end

function cluster_r2(r1s::Array{String,1},r2_seqs::Dict{String,String};id=0.97,t=12)
    tmpdir = mktempdir("/dev/shm"; prefix="vsearch_", cleanup=false)
    f_names = "$tmpdir/input.fasta"
    r_names = "$tmpdir/clustered.fasta"
    best = "$tmpdir/best.fasta"
    fo = open(f_names,"w")
    for n in r1s
        seq = r2_seqs[n]
        println(fo,">$n\n$seq\n")
    end
    close(fo)
    cluster = run(`vsearch --quiet --cluster_size $f_names --sizeout --threads $t --centroids $r_names --id $id`)
    cluster = run(`vsearch --quiet --sortbysize $r_names  --output $best --topn 1`)
    out = ""
    FastaReader(best) do fr
        for (n, seq) in fr
            out = String(seq)
        end
    end
    rm(tmpdir, recursive=true)
    return(out)
end




"Reads in fasta file in a dictionary name - sequencd"
function parse_fastq(f::String)::Dict{String,FASTX.FASTQ.Record}
    out = Dict{String, FASTX.FASTQ.Record}()
    ct = 0 
    for record in FASTQ.Reader(GzipDecompressorStream(open(f)))
        ct +=1
        id = FASTX.FASTQ.identifier(record)
        out[id]=record 
        if mod(ct,10) == 10
            println(stderr,"Parsed $ct sequences\r")
        end
    end
    println(stderr, "Parsed $ct sequences from $f")
    return(out)
end

function parse_mapping(mapfile::String)
    pairs=Array{Tuple{String,String},1}()
    r2c=Dict{String,String}() #read to cluster representative
    c2r=Dict{String,Array{String,1}}() #cluster rep to reads map
    for l in readlines(mapfile)
        parts = split(l)
        push!(pairs,(parts[1],parts[2]))
        r2c[String(parts[1])]=String(parts[2])
    end
    sort!(pairs, lt=(x,y)->isless(x[2],y[2]))
    cur_r2 = ""
    len = length(pairs)
    ct = 0
    cl_ct = 0
    r1s = Array{String,1}()
    cluster_4_work = Array{Tuple{Array{String,1},String},1}()
    @showprogress 1 "Parsing mapping..."  for pair in pairs
        ct += 1
        if ct == 1
            cur_r2 = pair[2]
        end
        if pair[2] != cur_r2 && ct !=1
            cl_ct += 1
            push!(cluster_4_work,(r1s,cur_r2))
            r1s = Array{String,1}()
            cur_r2 = pair[2]
        end
        c2r[cur_r2] = r1s 
        push!(r1s,pair[1])
    end
    push!(cluster_4_work,(r1s,cur_r2))
    c2r[cur_r2] = r1s 
    return(cluster_4_work,r2c,c2r)
end

function intersect_R1_R2_clusters(gjcR1,gjcR2)
    cluster_sets_r1, r2c_r1, c2r_r1 = parse_mapping(gjcR1)
    cluster_sets_r2, r2c_r2, c2r_r2 = parse_mapping(gjcR2)
    intersected_clusters = Array{Array{String,1},1}()
    @showprogress 1 "Geting R1 and R2 clustering intersections..." for r1c in keys(c2r_r1)
        #collect matching clusters i r2
        r2cs = Array{String,1}()
        for r1csR in c2r_r1[r1c] 
            push!(r2cs,r2c_r2[r1csR])
        end 
        r2cs = unique(r2cs)
        for r2_centroid in r2cs 
            r2_cluster = c2r_r2[r2_centroid]
            r2_cluster_intersect = intersect(r2_cluster,c2r_r1[r1c])
            push!(intersected_clusters,r2_cluster_intersect)
        end 
    end 
    println(stderr, "R1 clusters: ",length(c2r_r1)," R2 clusters: ",length(c2r_r2)," Combined clustering: ",length(intersected_clusters) )
    return(intersected_clusters)
end 

function julia_main()::Cint
    parser = ArgumentParser(description = "Rereplicate a clustered fasta sequence - for repairing...")

    add_argument!(parser, "-1", "--input-fastq1", type = String, help = "Input forward (R1) fastq")
    add_argument!(parser, "-2", "--input-fastq2", type = String, help = "Input reverese (R2) fastq")
    add_argument!(parser, "-a", "--clustering-r1", type = String, help = "A gjc of R1 file ")
    add_argument!(parser, "-b", "--clustering-r2", type = String, help = "A gjc of R2 file ")
    add_argument!(parser, "-o", "--output-fasta", type = String, help = "Output file")
    args = parse_args(parser)
    r2_fastq = parse_fastq(args.input_fastq2)    
    combined_clusters =  intersect_R1_R2_clusters(args.clustering_r1,args.clustering_r2)
    new_centroids = run_clustering(r2_fastq, combined_clusters)
    println(stderr, "Number of new centroids: ",length(new_centroids))
    exit()
    seqs = parse_fastq(args.input_fastq)
    expected_clustering_filei1 = args.input_fasta*".gjc"  

    if isfile(expected_clustering_file)
        println("Reading clustering from $expected_clustering_file")
    else 
        error("Expected clustering file $expected_clustering_file not found")
    end 
    write_out_fastq(seqs,mapping,args.output_fasta)
    return 0
end

julia_main()

