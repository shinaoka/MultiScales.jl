using Test
using MultiScales
using ITensors

@testset "imaginarytime.jl" begin
    @testset "decompose" begin
        β = 10.0
        ω = 0.1
        nbit = 10
        nτ = 2^nbit

        gtau(τ) = - exp(-ω * τ) / (1 + exp(-ω * β))
        @assert gtau(0.0) + gtau(β) ≈ -1

        τs = collect(LinRange(0.0, β, nτ + 1))[1:end-1]
        gtau_smpl = Vector{ComplexF64}(gtau.(τs))

        sites = siteinds("Qubit", nbit)
        gtau_mps = MultiScales.decompose(gtau_smpl, sites; cutoff = 1e-20)

        #@show gtau_mps

        gtau_smpl_reconst = vec(Array(reduce(*, gtau_mps), sites...))

        @test gtau_smpl_reconst ≈ gtau_smpl
    end
    #@test false
end
