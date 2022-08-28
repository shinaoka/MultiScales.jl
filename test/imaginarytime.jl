using Test
using MultiScales
import ITensors: siteinds
import SparseIR: Fermionic, Bosonic, FermionicFreq, valueim

function _test_data_imaginarytime(nbit, β)
    ω = 0.5
    N = 2^nbit
    halfN = 2^(nbit-1)

    # Tau
    gtau(τ) = - exp(-ω * τ) / (1 + exp(-ω * β))
    @assert gtau(0.0) + gtau(β) ≈ -1
    τs = collect(LinRange(0.0, β, N + 1))[1:end-1]
    gtau_smpl = Vector{ComplexF64}(gtau.(τs))

    # Matsubra
    giv(v::FermionicFreq) = 1/(valueim(v, β) - ω)
    vs = FermionicFreq.(2 .* collect(-halfN:halfN-1) .+ 1)
    giv_smpl = giv.(vs)

    return gtau_smpl, giv_smpl
end

@testset "imaginarytime.jl" begin
    @testset "decompose" begin
        β = 2.0
        nbit = 10
        nτ = 2^nbit

        gtau_smpl, giv_smpl = _test_data_imaginarytime(nbit, β)

        sites = siteinds("Qubit", nbit)
        gtau_mps = MultiScales.decompose(gtau_smpl, sites; cutoff = 1e-20)

        #@show gtau_mps

        gtau_smpl_reconst = vec(Array(reduce(*, gtau_mps), reverse(sites)...))

        @test gtau_smpl_reconst ≈ gtau_smpl
    end

    @testset "ImaginaryTimeFT" begin
        ITensors.set_warn_order(100)
        β = 1.5
        nbit = 6
        nτ = 2^nbit

        gtau_smpl, giv_smpl = _test_data_imaginarytime(nbit, β)

        sites = [Index(2, "Qubit,τ=$t,iω=$(nbit+1-t)") for t in 1:nbit]
        gtau_mps = MultiScales.decompose(gtau_smpl, sites; cutoff = 1e-20)
        ft = MultiScales.ImaginaryTimeFT(MultiScales.FTCore(sites))
        giv_mps = MultiScales.to_wn(Fermionic(), ft, gtau_mps, β; cutoff = 1e-20)

        # w_Q, ..., w_1
        giv =  vec(Array(reduce(*, giv_mps), sites...))

        #open("test.txt", "w") do file
            #for i in 1:nτ
                #println(file, i, " ", real(giv[i]), " ", imag(giv[i]), " ", real(giv_smpl[i]), " ", imag(giv_smpl[i]))
            #end
        #end
        ##@show giv_mps

        @test maximum(abs, giv - giv_smpl) < 2e-2
        #@test false
    end

    #@test false
end
