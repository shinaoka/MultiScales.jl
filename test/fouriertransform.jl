using Test
using MultiScales
using ITensors

@testset "fouriertransform.jl" begin
    @testset "_qt" for sign in [-1, 1]
        nbit = 4
        N = 2^nbit
        
        sites = siteinds("Qubit", nbit)
        M = MultiScales._qft(sites; sign=sign)

        # Return the bit of an integer `i` at the position `pos` (`pos=1` is the least significant digit).
        bitat(i, pos) = ((i & 1<<(pos-1))>>(pos-1))
        @assert bitat(2, 1) == 0
        @assert bitat(2, 2) == 1

        # Input function `f(x)` is 1 only at xin otherwise 0.
        xin = 3 # 1-start
        @assert xin <= N
        tmp = collect(string(bitat(xin-1, pos)) for pos in nbit:-1:1)
        mpsf = MPS(sites, tmp)

        # Physical indices from lef to right: t_Q, ..., t_1
        mpsg = reduce(*, noprime(contract(M, mpsf)))

        # Values of output function
        outfunc = vec(Array(mpsg, sites))

        @test outfunc ≈ [exp(sign * im * 2π * (t-1) * (xin-1)/N)/sqrt(N) for t in 1:N]
    end
end
