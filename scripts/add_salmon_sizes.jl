using ArgParse2
using CSV
using DataFrames
using DataStructures
using BioSequences
using CodecZlib
using BioAlignments
using DataFramesMeta

"Reads in fasta file in a dictionary name - sequencd"
function parse_fasta(f::String)::Dict{String,String}
    out = Dict{String,String}()
    ct = 0
    for record in FASTA.Reader(open(f))
        ct +=1
        id = BioSequences.FASTA.identifier(record)
        seq = BioSequences.FASTA.sequence(record)
        out[id]=String(seq)
        if mod(ct,10) == 0
            println(stderr,"Parsed $ct sequences\r")
        end
    end
    println(stderr, "Parsed $ct sequences from $f")
    return(out)
end

function get_salmon_sizes(infile::String,duplicates::String)::Dict{String,String}
    salmon_data = DataFrame(CSV.read(infile))
    salmon_duplicates = DataFrame(CSV.read(duplicates))
    mintpm = minimum(salmon_data.TPM) + 1
    salmon_data.salmon_size = Int64.(round.(salmon_data.NumReads)) .+ 1
    sort!(salmon_data,:TPM)
    out = Dict{String,String}()
    for r in eachrow(salmon_data)
	name = split(r.Name,";")[1]
	ns = r.salmon_size
	new_name = "$name;size=$ns"
	out[r.Name] = new_name
    end
    for r in eachrow(salmon_duplicates)
	id_dup = split(r.DuplicateRef,";")[1]
	id_ref = split(r.RetainedRef,";")[1]
	out[r.DuplicateRef] = replace(out[r.RetainedRef],id_ref => id_dup)
    end
    return(out)
	
end

function main()
    parser = ArgumentParser(prog = "add_salmon_sizes.jl",
                        description = "Add salmon TPM estimates as pseudocounts on vsearch format to contigs",
                        epilog = "sizes..ajusted")
    add_argument!(parser, "--infasta", "-i", help = "input fasta ",type = String)
    add_argument!(parser, "--outfasta", "-o", help = "output fasta ",type = String)
    add_argument!(parser, "--salmon-counts","-s", help = "salmon counts",type = String)
    add_argument!(parser, "--salmon-duplicates","-d", help = "salmon duplicates",type = String)
    args = parse_args(parser)

    new_names = get_salmon_sizes(args.salmon_counts,args.salmon_duplicates)
    outf = open(args.outfasta,"w")
    for record in FASTA.Reader(open(args.infasta))
        id = BioSequences.FASTA.identifier(record)
        seq = BioSequences.FASTA.sequence(record)
	outid = new_names[id]
	print(outf,">$outid\n$seq\n")
    end
    close(outf)
end

main()
