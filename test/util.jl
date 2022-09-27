using Test
import MultiScales
using ITensors

@testset "util.jl" begin
    @testset "replace_mpo_siteinds!" begin
        nbit = 3
        sites = siteinds("Qubit", nbit)
        M = MPO(sites, ["Y" for n in 1:nbit])
        #@show sites

        sites2 = [Index(2, "n=$n") for n in 1:nbit]
        MultiScales.replace_mpo_siteinds!(M, sites, sites2)

        # Tests fail weird way if @test is used.
        @assert all([!hasind(M[n], sites[n]) for n in 1:nbit])
        @assert all([!hasind(M[n], sites[n]') for n in 1:nbit])
        @assert all([hasind(M[n], sites2[n]) for n in 1:nbit])
        @assert all([hasind(M[n], sites2[n]') for n in 1:nbit])
        true
    end

    @testset "splitsiteind" begin
        nbit = 6
        sites = siteinds("Qubit", nbit)
        csites = [Index(4, "csite=$s") for s in 1:nbit÷2]
        D = 3
        mps = randomMPS(csites; linkdims=D)
        mps_split = MultiScales.splitsiteind(mps, sites)
        @test vec(Array(reduce(*, mps_split), sites)) ≈ vec(Array(reduce(*, mps), csites))

        mps_reconst = MultiScales.combinesiteinds(mps_split, csites)
        @test vec(Array(reduce(*, mps_reconst), csites)) ≈ vec(Array(reduce(*, mps), csites))
    end
end