using ArgParse2
using BioSequences
using CodecZlib
using FastaIO
using ProgressMeter
using Base.Threads



function write_out(centroids::Array{String,1},mapping::Array{Tuple{Array{String,1},String},1},outf::String)
    fo = open(outf,"w")
    println("Riting out to $outf  ...")
    @showprogress 1 "Writing out to $outf ..." for i in 1:length(centroids)
        seq = centroids[i]
        name = mapping[i][2]
        print(fo,">$name\n$seq\n")
    end
    close(fo)
end

function write_out_fastq(seqs::Dict{String,BioSequences.FASTQ.Record}, mapping::Array{Tuple{Array{String,1},String},1},outf::String)
    fo = FASTQ.Writer(GzipCompressorStream(open(outf,"w")))
    println("Riting out to $outf  ...")
    @showprogress 1 "Writing out to $outf ..." for pair in mapping 
        centroid = seqs[pair[2]]
        for seq_id in pair[1]
            cluster_mem =  FASTQ.Record(seq_id, BioSequences.FASTQ.sequence(centroid), BioSequences.FASTQ.quality(centroid))
            write(fo,cluster_mem)
        end 
        #seq = centroids[i]
        #name = mapping[i][2]
        #print(fo,">$name\n$seq\n")
    end
    close(fo)
end



function run_clustering(seqs::Dict{String,String}, mapping::Array{Tuple{Array{String,1},String},1})
    #mapping = mapping[1:10000]
    n = length(mapping)
    p = Progress(n)
    println(stderr,"Starting $n clustering runs...")
    centroids = Array{String,1}(undef, n)
    Threads.@threads for i in 1:n
                r1s = mapping[i][1]
                name = mapping[i][2]
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
function parse_fastq(f::String)::Dict{String,BioSequences.FASTQ.Record}
    out = Dict{String, BioSequences.FASTQ.Record}()
    ct = 0 
    for record in FASTQ.Reader(GzipDecompressorStream(open(f)))
        ct +=1
        id = BioSequences.FASTQ.identifier(record)
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
    for l in readlines(mapfile)
        parts = split(l)
        push!(pairs,(parts[1],parts[2]))
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
        push!(r1s,pair[1])
    end
    push!(cluster_4_work,(r1s,cur_r2))
    return(cluster_4_work)
end

function julia_main()::Cint
    parser = ArgumentParser(description = "Rereplicate a clustered fasta sequence - for repairing...")

    add_argument!(parser, "-f", "--input-fasta", type = String, help = "Input fastq")
    add_argument!(parser, "-q", "--input-fastq", type = String, help = "Fatq matching fasta after clustering clustered Must have also a gjc file ")
    add_argument!(parser, "-o", "--output-fasta", type = String, help = "Output file")
    args = parse_args(parser)
    seqs = parse_fastq(args.input_fastq)
    expected_clustering_file = args.input_fasta*".gjc"  
    mapping = parse_mapping(expected_clustering_file)

    if isfile(expected_clustering_file)
        println("Reading clustering from $expected_clustering_file")
    else 
        error("Expected clustering file $expected_clustering_file not found")
    end 
    write_out_fastq(seqs,mapping,args.output_fasta)
    return 0
end

julia_main()
