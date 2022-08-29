using Test
import MultiScales
import ITensors: siteinds, Index

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
end