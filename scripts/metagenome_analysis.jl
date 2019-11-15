#!/usr/bin/env julia
using ArgParse
using CSV
using DataFrames
using Query 

function analyse(files::Array{String,1}, cut_off::Float64, output_stem::String)
  data = DataFrame(Taxon=Array{String,1}(),Fraction_of_reads=Array{Float64,1}(),Sample=Array{String,1}())
  DI = DataFrame(Sample=Array{String,1}(),Shannon_DI=Array{Float64,1}(),N_1proc=Array{Int64,1}(),N_5proc=Array{Int64,1}())
  DI_part = 0.0
  ct_1 = 0
  ct_5 = 0
  ct_all = 0

  for f in files
    data_p = CSV.read(f,delim="\t")
    sample = split(basename(f),'.')[1]
    ct_1 = 0
    ct_5 = 0
    ct_all = 0
    DI_part = 0
    for r in eachrow(data_p)
      p = r[:fraction_total_reads]
      ct_all += 1
      if p >= 0.01
        ct_1 += 1
        if p >= 0.05
          ct_5 += 1
        end
      end
      if p > 0
        DI_part += p*log(p)
      end
      push!(data,[r[:name],p,sample])
    end 
    push!(DI,[sample,-1*DI_part,ct_1,ct_5])
  end
  #filter out low counts
  data = @from i in data  begin
            @where i.Fraction_of_reads >= cut_off 
            @select {i.Taxon, i.Fraction_of_reads,i.Sample}
            @collect DataFrame
       end
  sort!(data, [:Sample, :Fraction_of_reads])
  CSV.write(output_stem*"_filtered_fractions.csv",data)
  CSV.write(output_stem*"_additional_info.csv",DI)
end


function main(args)

    #= Command-line option parser =#
    arg_parse_settings = ArgParseSettings(description="Program calculates pairs of reads with not matching names inserts")
    @add_arg_table arg_parse_settings begin
        "--abundances", "-r"
        nargs = '+'
        arg_type = String
        help = "BRACKEN outputs"
        required = true
        "--cut-off", "-c"
        arg_type = Float64 
        help = "Cut off of abundances (part of unit)"
        required = true
        "--output", "-o"
          arg_type = String
          help = "Output stem"
          required = true
     end

    parsed_args = parse_args(arg_parse_settings) #= In order to use variables, =#
                                                 #= use parsed_args["foo_bar"] =#
    #= Main code =#
    analyse(parsed_args["abundances"], parsed_args["cut-off"],parsed_args["output"])
end

main(ARGS)
