
using ArgParse2
using CSV
using DataFrames
using DataStructures
using BioSequences
using CodecZlib
using BioAlignments
using DataFramesMeta

"Reads in fasta file in a dictionary name - sequencd"
function parse_fasta(f::String, o::String)
    ct = 0
    fo = open(o,"w")
    for record in FASTA.Reader(open(f))
        ct +=1
        id = BioSequences.FASTA.identifier(record)
        idgene = split(id,"_")[1]
        seq = BioSequences.FASTA.sequence(record)
        seq_len = length(seq)
        gtfl = """$id\tASSEMBLY\texon\t1\t$seq_len\t0.000000\t+\t.\tgene_id "$idgene"; transcript_id "$id"; transcript_biotype "Contig"; gene_biotype "Contig";"""
        println(fo, gtfl)
    end
    close(fo)
end

function main()
    parser = ArgumentParser(prog = "analyse_alignment_on_reference.jl",
                        description = "Evaluate how many 16s copies we are able to identify",
                        epilog = "Removing potential chimeras and repair snps")
    add_argument!(parser, "--infasta", "-i", help = "input sam ",type = String)
    add_argument!(parser, "--outgtf", "-o", help = "Sumary statistics ",type = String)
    args = parse_args(parser)

    parse_fasta(args.infasta, args.outgtf)
end

main()
