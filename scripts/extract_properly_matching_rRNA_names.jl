using ArgParse
using CSV
using DataFrames
using DataStructures
using BioSequences
using CodecZlib


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
            help = "Input file"
            required = true 
        "--log", "-l"
            help = "Log file to write reads on target"
            arg_type = String
            required = true 
        "--names", "-n"
            help = "Write out names of reads "
            arg_type = String
            required = true 
        "--fastqin", "-q"
            help = "Fasq to filter out reads"
            arg_type = String
            required = true 
        "--fastqout", "-o"
            help = "Fasqout to filter out reads"
            arg_type = String
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

function main()
    parsed_args = parse_commandline()
    logf = open(parsed_args["log"],"w")
    namesf = open(parsed_args["names"],"w")
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
    println(stderr,"Filtering out matches in template between $t_star:$t_end")
    chosen_reads = Set{String}() 
    # read the hmm output 
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
                    ts, te = sort([ts,te])
                    rs, re = sort([rs,re])
                    #println("$r_star $r_end $t_star $t_end   real  $rs $re  $ts $te ")
                    if ( (re-rs+1) >= minl ) && (ts >= t_star) && (ts <= t_end) && (te >= t_star) && (te <= t_end) && (rs >= r_star) && (rs <= r_end) && (re >= r_star) && (re <= r_end)
                        println(namesf,rname)
                       push!(chosen_reads,rname)
                    end
                end
                line += 1
            end
    end
    close(namesf)
    ct_t = length(chosen_reads)
    ct_p = 0
    ct_all = 0
    #read in fastq and filter
    fastq_file = parsed_args["fastqin"] 
    fastq_fileo = parsed_args["fastqout"] 
    reader = FASTQ.Reader(GzipDecompressorStream(open(fastq_file)))
    writer = FASTQ.Writer(GzipCompressorStream(open(fastq_fileo,"w")))
    println("Choosing reads with proper 16S fragments in file: $fastq_file")
    record = FASTQ.Record()
    while !eof(reader)
        read!(reader, record)
        ct_all += 1
        #println(record)
        id = FASTQ.identifier(record)
        if id in chosen_reads 
            ct_p += 1
            write(writer,record) 
            delete!(chosen_reads,id)
        end
        ## Do something.
    end
    if ct_p != ct_t
        println(stderr, "Something wrong in files: $ct_t in csv and $ct_p in the fastq filr")
    end
    frac = ct_p/ct_all
    println(stderr,"On 'target' $frac reads")
    println(logf,"$frac")
    close(writer)



    
end

main()
