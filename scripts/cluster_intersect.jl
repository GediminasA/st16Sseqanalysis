using BioSequences
using CodecZlib
using FastaIO
using FASTX 
using Base.Threads
using ProgressMeter
using Suppressor
push!(LOAD_PATH, "scripts/julia_modules/ArgParse2.jl/src")
println()
using ArgParse2



function write_out(centroids::Array{String,1},mapping::Array{Tuple{Array{String,1},String},1},outf::String)
    fo = open(outf,"w")
    println("Writing out to $outf  ...")
    @showprogress 1 "Writing out to $outf ..." for i in 1:length(centroids)
    #for i in 1:length(centroids)
        seq = centroids[i]
        name = mapping[i][2]
        print(fo,">$name\n$seq\n")
    end
    close(fo)
end

function write_out_fastq(seqs1::Dict{String,FASTX.FASTQ.Record},
                         seqs2::Dict{String,FASTX.FASTQ.Record},
                         mapping::Array{String,1},
                         outf1::String,
                         outf2::String,
                        )
    fo1 = FASTQ.Writer(GzipCompressorStream(open(outf1,"w")))
    fo2 = FASTQ.Writer(GzipCompressorStream(open(outf2,"w")))
    println("Writing out to $outf1  ...")
    println("Writing out to $outf2  ...")
    @showprogress 1 "Writing out to $outf1 $outf2 ..." for seq_id in mapping 
        centroid1 = seqs1[seq_id]
        centroid2 = seqs2[seq_id]
        cluster_mem1 =  FASTQ.Record(seq_id, FASTX.FASTQ.sequence(centroid1), FASTX.FASTQ.quality(centroid1))
        cluster_mem2 =  FASTQ.Record(seq_id, FASTX.FASTQ.sequence(centroid2), FASTX.FASTQ.quality(centroid2))
        write(fo1,cluster_mem1)
        write(fo2,cluster_mem2)
    end
    close(fo1)
    close(fo2)
end

function run_clustering(seqs::Dict{String,FASTX.FASTQ.Record}, mapping::Array{Array{String,1},1};swarm=false)
    seqs2 = Dict{String,String}()
    for k in keys(seqs)
        seqs2[k] = FASTX.FASTQ.sequence(seqs[k])
    end 
    run_clustering(seqs2,mapping,swarm=swarm)
end
 
function run_clustering(seqs::Dict{String,String}, mapping::Array{Array{String,1},1};swarm=false)
    #mapping = mapping[1:10000]
    n = length(mapping)
    p = Progress(n)
    println(stderr,"Starting $n clustering runs...")
    centroids = Array{String,1}(undef, n)
    Threads.@threads for i in 1:n
                r1s = mapping[i]
                if length(r1s) == 1
                    centroids[i] = r1s[1]
                else
                    centroids[i] = cluster_r2(r1s, seqs,t=3,swarm=swarm)
                end
                next!(p)
    end
    return(centroids)
end

function cluster_r2(r1s::Array{String,1},r2_seqs::Dict{String,String};id=0.94,t=12,swarm=false)
    tmpdir = mktempdir("/dev/shm"; prefix="vsearch_", cleanup=false)
    f_names = "$tmpdir/input.fasta"
    f_names_de = "$tmpdir/input_woident.fasta"
    f_names_de_won = "$tmpdir/input_woident_won.fasta"
    r_names = "$tmpdir/clustered.fasta"
    best = "$tmpdir/best.fasta"
    fo = open(f_names,"w")
    for n in r1s
        seq = r2_seqs[n]
        println(fo,">$n\n$seq\n")
    end
    close(fo)
    cluster = run(`vsearch --quiet --derep_fulllength $f_names --sizeout  --output $f_names_de`)
    if swarm
        try 
            @suppress begin
                cluster = run(`bbduk.sh  maxns=0  in=$f_names_de out=$f_names_de_won`)
                cluster = run(`swarm -z  $f_names_de_won -t $t -w $r_names -d 2 `)
            end 
        catch 
            cluster = run(`vsearch --quiet --sizein --cluster_size $f_names_de --sizeout --threads $t --centroids $r_names --id $id`)
        end 
    else 
        cluster = run(`vsearch --quiet --sizein --cluster_size $f_names_de --sizeout --threads $t --centroids $r_names --id $id`)
    end 
    cluster = run(`vsearch --quiet --sortbysize $r_names  --output $best --topn 1`)
    out = ""
    FastaReader(best) do fr
        for (n, seq) in fr
            out = split(String(n),";")[1]
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
    add_argument!(parser, "-o", "--output-fastq1", type = String, help = "Output file 4 R1")
    add_argument!(parser, "-p", "--output-fastq2", type = String, help = "Output file 4 R2")
    add_argument!(parser,"-s","--swarm", action = "store_true",help="Use swarm instead of vsearch for clustering")
    add_argument!(parser, "-m", "--min-cluster-size", type = Int, default=1 , help = "Minimum cluster size")
    args = parse_args(parser)
    r1_fastq = parse_fastq(args.input_fastq1)    
    r2_fastq = parse_fastq(args.input_fastq2)    
    combined_clusters_ini =  intersect_R1_R2_clusters(args.clustering_r1,args.clustering_r2)
    #filter based on size
    combined_clusters = Array{Array{String,1},1}()
    for c in combined_clusters_ini  
        if length(c) >= args.min_cluster_size 
            push!(combined_clusters,c)
        end 
    end 
    println(stderr, "Number of cluster after size filtering ",length(combined_clusters))
    new_centroids = run_clustering(r2_fastq, combined_clusters, swarm=args.swarm)
    println(stderr, "Number of new centroids: ",length(new_centroids))
    write_out_fastq(r1_fastq, r2_fastq, new_centroids, args.output_fastq1, args.output_fastq2)        
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

