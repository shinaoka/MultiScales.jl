using Test
import MultiScales
using ITensors

function _tomat(a)
    sites = siteinds(a)
    N = length(sites)
    halfN = N ÷ 2
    sites_ = [sites[1:2:N]..., sites[2:2:N]...]
    return reshape(Array(reduce(*, a), sites_), 2^halfN, 2^halfN)
end

@testset "matmul.jl" begin
    @testset "matmul" begin
        nbit = 6
        sites = siteinds("Qubit", nbit)
        csites = [Index(4, "csite=$s") for s in 1:nbit÷2]

        D = 3
        a = randomMPS(sites; linkdims=D)
        b = randomMPS(sites; linkdims=D)

        ab = MultiScales.matmul(a, b)

        abmat = _tomat(a) * _tomat(b)

        @test _tomat(ab) ≈ abmat
    end
    #@test false
end