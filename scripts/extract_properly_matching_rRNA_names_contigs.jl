using ArgParse
using CSV
using DataFrames
using BioSequences

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--16sr", "-r"
            help = "16S expected region in a read, <start>:<end> "
            arg_type = String
            required = true 
        "--16st", "-t"
            help = "16S expected region in a hmm/alignment, <start>:<end> "
            arg_type = String
            required = true 
        "--min", "-m"
            help = "minimum match length "
            arg_type = Int 
            required = true 
        "--input","-i"
            help = "Input file hmmscan"
            required = true 
        "--fasta","-f"
            help = "Input file hmmscan"
            required = true 
        "--matching-fasta","-o"
            help = "Input file hmmscan"
            required = true 
        "--matching-motifs","-g"
            help = "Input file hmmscan"
            required = true 

        # "--flag1"
        #     help = "an option without argument, i.e. a flag"
        #     action = :store_true
        # "arg1"
        #     help = "a positional argument"
        #     required = true
    end

    return parse_args(s)
end

function read_fasta(inf::String)::Dict{String,String}
    reader = open(FASTA.Reader, inf)
    record = FASTA.Record()
    out = Dict{String,String}()
    while !eof(reader)
        read!(reader, record)
        seq::String = FASTA.sequence(record)   
        name = FASTA.identifier(record)
        out[name] = seq 
    end 
    return(out)
end

function main()
    parsed_args = parse_commandline()
    #output
    out_s = open(parsed_args["matching-fasta"],"w")
    out_m = open(parsed_args["matching-motifs"],"w")
    #read fasta

    seqs = read_fasta(parsed_args["fasta"])
    # get read limit
    r_data= split(parsed_args["16sr"],":")
    r_star = parse(Int64,r_data[1])
    r_end = parse(Int64,r_data[2])
    println(stderr,"Filtering out matches in reads between $r_star:$r_end")
    
    # get red limit
    t_data= split(parsed_args["16st"],":")
    t_star = parse(Int64,t_data[1])
    t_end = parse(Int64,t_data[2])
    t_star, t_end = sort([t_star,t_end])
    r_star, r_end = sort([r_star,r_end])

    minl = parsed_args["min"]
    outputed = Array{String,1}()
    println(stderr,"Filtering out matches in template between $t_star:$t_end")
    open(parsed_args["input"]) do f
            line = 1
            while !eof(f)
                x = readline(f)
                if x[1] != '#'
                    parts = split(x)
                    rname = parts[1]
                    ts = parse(Int64,parts[5])
                    te = parse(Int64,parts[6])
                    rs = parse(Int64,parts[9])
                    re = parse(Int64,parts[10])
                    rlen = parse(Int64,parts[11])
                    strand = parts[12] 
                    ts, te = sort([ts,te])
                    rs, re = sort([rs,re])
                    if strand == "+"
                        rs, re = sort(([rs,re] .- rlen .- 1) * -1 )
                        ss = reverse_complement(BioSequence{DNAAlphabet{4}}(seqs[rname]))
                        seqs[rname] = string(ss)
                    end
                    if ( (re-rs+1) >= minl ) && (ts >= t_star) && (ts <= t_end) && (te >= t_star) && (te <= t_end) && (rs >= r_star) && (rs <= r_end) && (re >= r_star) && (re <= r_end) && !(rname in outputed)
                        push!(outputed,rname)
                        motiv_seq = seqs[rname][rs:re]
                        seq = seqs[rname]
                        println(out_m,">$rname\n$motiv_seq\n")
                        println(out_s,">$rname\n$seq\n")

                        println(stdout,rname)
                    end
                end
                line += 1
            end
    end
    close(out_s)
    close(out_m)
    # read the hmm output 
    
end

main()
