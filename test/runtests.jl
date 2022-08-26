using Test
using MultiScales

@testset "MultiScales.jl" begin
    include("arithmetic.jl")
    include("fouriertransform.jl")
end
