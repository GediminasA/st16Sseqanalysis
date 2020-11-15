push!(LOAD_PATH, dirname(Base.PROGRAM_FILE))
using sanitise
using ArgParse2
function julia_main()::Cint
    parser = ArgumentParser(prog = "tidy_assembly.jl",
                            description = "Analyse assembly of contigs and prepares it for clustering",
                            epilog = "Removing potential chimeras and repair snps")

    add_argument!(parser, "-i", "--input-fasta", type = String, help = "Input file")
    add_argument!(parser, "-o", "--output-fasta", type = String, help = "Output file")
    add_argument!(parser, "-w", "--windows-size", type = Int64, help = "Smoothing windows size")
    add_argument!(parser, "-s", "--windows-step", type = Int64, help = "Smoothing windows step")
    add_argument!(parser, "-c", "--smoothing-cutoff", type = Float64, help = "Smoothing windows cutoff")
    args = parse_args(parser)
    aln = sanitise.Alignment(args.input_fasta)
    vsearh_cluster!(aln,args.windows_size,args.windows_step,args.smoothing_cutoff)
    write_to_fasta(aln,args.output_fasta)

    return 0
end

julia_main()
#mycoolfunction()
