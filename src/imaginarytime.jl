@doc """
For imaginary-time/-frequency domains
"""
struct ImaginaryTimeFT <: AbstractFT
    ftcore::FTCore

    function ImaginaryTimeFT(ftcore::FTCore)
        new(ftcore)
    end
end


function to_wn(::Type{Fermionic}, ft::ImaginaryTimeFT, gtau::MPS; kwargs...)
    length(gtau) == nbit(ft) || error("Length mismatch")
    gtau = noprime(copy(gtau))

    sites = [n == 1 ? ind(gtau[n], 1) : ind(gtau[n], 2) for n in eachindex(gtau)]

    # Apply phase shift to each Qubit
    θ = π * ((-N+1)/N)
    for i in 1:nbit(ft)
        gtau[i] *= noprime(gtau[i] * op("Phase", sites[i]; ϕ= θ * 2^(nbit-i) ))
    end

    # FFT
    return apply(ft.ftcore.forward, kwargs...)
end


function decompose(gtau_smpl::Vector{ComplexF64}, sites; kwargs...)
    nbit = length(sites)
    length(gtau_smpl) == 2^nbit || error("Length mismatch")

    # (g_1, g_2, ...)
    gtau_smpl = reshape(gtau_smpl, repeat([2,], nbit)...)
    permutedims!(gtau_smpl, gtau_smpl, reverse(collect(1:nbit)))

    return MPS(gtau_smpl, sites; kwargs...)
end