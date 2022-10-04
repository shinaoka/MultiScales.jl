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

        abmat = _tomat(a) * _tomat(b)

        @test _tomat(MultiScales.matmul(a, b)) ≈ abmat
        @test _tomat(MultiScales.matmul_naive(a, b)) ≈ abmat
    end

    @testset "matmul_thru_mpo" begin
        nbit = 6
        sites = siteinds("Qubit", nbit)
        csites = [Index(4, "csite=$s") for s in 1:nbit÷2]

        D = 2
        a = randomMPS(sites; linkdims=D)
        b = randomMPS(sites; linkdims=D)
        mpo_a = MultiScales.tompo_matmul(a, csites)

        b_ = MultiScales.combinesiteinds(b, csites)
        ab = apply(mpo_a, b_)
        ab = MultiScales.splitsiteind(ab, sites)

        abmat = _tomat(a) * _tomat(b)
        @test _tomat(ab) ≈ abmat
    end
end