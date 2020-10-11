push!(LOAD_PATH, "../../")
using sanitise
module Testsanitise

using Test



# include("AlignmentMatrix.jl")
# exit()
 @testset "AlignmentMatrix" begin
     include("AlignmentMatrix.jl")
 end

# @testset "Utilities" begin
#     include("utilities.jl")
# end



end #module TestOTDDNprimer
