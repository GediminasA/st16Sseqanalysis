using ArgParse
using CSV
using DataFrames
using DataStructures
using BioSequences
using CodecZlib
using BioAlignments
using BGZFStreams

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--inbam", "-i"
            help = "input sam "
            arg_type = String
            required = true 
        "--outbam", "-o"
            help = "Filtered sam "
            arg_type = String
            required = true 
        "--mintlen", "-l"
            help = "minimum TLEN  length "
            arg_type = Int 
            required = true 
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println(parsed_args["inbam"])
    reader = open(BAM.Reader,parsed_args["inbam"],index = parsed_args["inbam"]*".bai")
    h = BAM.header(reader)
    writer = BAM.Writer(BGZFStream(open(parsed_args["outbam"], "w"), "w"),h)
    for record in reader
        if BAM.ismapped(record) #&& SAM.isprimary(record)
            aln = BioAlignments.BAM.alignment(record)
            s = first(aln.anchors).refpos
            e = last(aln.anchors).refpos
            l = abs(s-e)+1
            if l >= parsed_args["mintlen"]
                write(writer,record)
            end
        end 
    end
    close(reader)
    close(writer)
    
end

main()
