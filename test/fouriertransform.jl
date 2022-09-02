using Test
using MultiScales
using ITensors

@testset "fouriertransform.jl" begin
    #===
    @testset "_qt_tensor" begin
        sign = 1

        # (nbit, layer, pos, sign)
        m =  MultiScales._tensor(1, 1, 1, sign)
        @test size(m) == (1, 2, 2)
        @test reshape(m, (2, 2)) ≈ [1 1; 1 -1]

        # (nbit=2, layer, pos, sign)
        m =  MultiScales._tensor(2, 1, 1, sign)
        @test false
    end
    ===#

    #@testset "_qt" for sign in [-1, 1]
    #@testset "_qt" for sign in [-1], nbit in [1], xin in 0:(2^nbit-1)
    @testset "_qt" for sign in [-1], nbit in [2], xin in [0]
        #nbit = 4
        N = 2^nbit
        
        sites = siteinds("Qubit", nbit)
        M = MultiScales._qft3(sites; sign=sign)
        #==
        for (i, m) in enumerate(M)
            println("")
            println("")
            println("")
            @show i, m
        end
        ==#

        # Return the bit of an integer `i` at the position `pos` (`pos=1` is the least significant digit).
        bitat(i, pos) = ((i & 1<<(pos-1))>>(pos-1))
        @assert bitat(2, 1) == 0
        @assert bitat(2, 2) == 1

        # Input function `f(x)` is 1 only at xin otherwise 0.
        #xin = 2 # 0-based
        #xin = 0 # 0-based
        @assert xin <= N - 1
        tmp = collect(string(bitat(xin, pos)) for pos in nbit:-1:1)
        @show tmp
        mpsf = MPS(sites, tmp)

        # Physical indices from left to right: x_Q, ..., x_1
        mpsg = reduce(*, noprime(contract(M, mpsf)))

        # Values of output function
        outfunc = vec(Array(mpsg, sites))

        @test outfunc ≈ [exp(sign * im * 2π * y * xin/N)/sqrt(N) for y in 0:(N-1)]
    end
    @test false
end
