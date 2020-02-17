using ArgParse
using CSV
using DataFrames
using DataStructures
using BioSequences
using CodecZlib
using BioAlignments
using BGZFStreams
using Statistics
using StatsBase

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--incsv", "-i"
            help = "input sam "
            arg_type = String
            required = true
        "--outstem", "-o"
            help = "stem for "
            arg_type = String
            required = true

    end

    return parse_args(s)
end


function main()

    parsed_args = parse_commandline()
    out_stem = parsed_args["outstem"]
    out_c = out_stem*"_first_letter_counts.csv"
    out_m = out_stem*"_insert_size_medians.csv"
    out_mall = out_stem*"_insert_size_medians_all.csv"
    out_h = out_stem*"_insert_sizes_histogram.csv"

    df = CSV.read(parsed_args["incsv"])

    # sarting leter count
    size_df = by(df,[:Species,:Starting_letter,:Sample], :Insert_size => length)
    size_df_unstacked_fl = unstack(size_df,:Starting_letter,:Insert_size_length)
    CSV.write(out_c, size_df_unstacked_fl)
    #length median
    size_df = by(df,[:Species,:Starting_letter,:Sample], :Insert_size => median)
    size_df_unstacked_median = unstack(size_df,:Starting_letter,:Insert_size_median)
    CSV.write(out_m, size_df_unstacked_median)
    #length median for all
    size_df = by(df,[:Sample], :Insert_size => median)
    CSV.write(out_mall, size_df)
    #
    #insert histogram
    df.id=1:size(df, 1)
    df_unstacked = unstack(df,:Sample,:Insert_size)
    samples = names(df_unstacked)[4:end]
    edges = 0.0:50.0:1500.0
    histogramdf = DataFrame()
    histogramdf[:Length] = collect(edges)[2:end]
    for s in samples
        d4h =  collect(skipmissing(df_unstacked[s]))
        res = StatsBase.fit(Histogram,d4h,edges)
        histogramdf[s] = res.weights ./ sum(res.weights)
    end
    println(histogramdf)
    CSV.write(out_h,histogramdf)


end

main()
