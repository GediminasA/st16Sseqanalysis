push!(LOAD_PATH, dirname(Base.PROGRAM_FILE)*"/julia_modules/")
using sanitise
using ArgParse2
function julia_main()::Cint
    parser = ArgumentParser(prog = "tidy_assembly.jl",
                            description = "Analyse assembly of contigs and prepares it for clustering",
                            epilog = "Removing potential chimeras and repair snps")

    add_argument!(parser, "-i", "--input-fasta", type = String, help = "Input file")
    args = parse_args(parser)
    aln = sanitise.Alignment(args.input_fasta)
    println(aln_conservation(aln,ignore_amounts=true))

    return 0
end

julia_main()
#mycoolfunction()
