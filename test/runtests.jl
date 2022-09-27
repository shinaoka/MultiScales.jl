using Test
using MultiScales

@testset "MultiScales.jl" begin
    include("util.jl")
    include("arithmetic.jl")
    include("fouriertransform.jl")
    include("imaginarytime.jl")
    include("matmul.jl")
end
