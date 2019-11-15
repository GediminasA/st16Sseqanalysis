using ArgParse
using CSV
using DataFrames

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
        "--max", "-m"
            help = "maximum match length "
            arg_type = Float64 
            required = true 
        "--input","-i"
            help = "Input file"
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

    maxl = parsed_args["max"]
    
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
                    rl = parse(Int64,parts[11])
                    ts, te = sort([ts,te])
                    rs, re = sort([rs,re])
                    limit = rl*maxl 
                    #println("$r_star $r_end $t_star $t_end   real  $rs $re  $ts $te ")
                    if ( (re-rs+1) >= limit) && (ts >= t_star) && (ts <= t_end) && (te >= t_star) && (te <= t_end) && (rs >= r_star) && (rs <= r_end) && (re >= r_star) && (re <= r_end)
                        println(stdout,rname)
                    end
                end
                line += 1
            end
    end
    # read the hmm output 
    
end

main()
