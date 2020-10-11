module sanitise
using Statistics
using StatsBase 
using FastaIO
using AutoHashEquals
using DataStructures
using DataFrames
export Alignment
export read_fasta
export write_to_fasta
export sub_alignment
export vsearh_cluster!
export get_max_cl_length
export char_array_entropy 
export aln_conservation
export read_jc_clusters
export mark_proper_size_clusters
export parse_dada_results 
export eval_recomb
include("AlignmentMatrix.jl")
include("utils.jl"	)
end
