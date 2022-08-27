abstract type AbstractFT end

struct FTCore
    forward::MPO

    function FTCore(nbit)
        sites = siteinds("Qubit", nbit)
        new(_qft(sites))
    end
end

nbit(ft::AbstractFT) = length(ft.ftcore.mpo)


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
    gtau = copy(gtau)
    # TODO: Check noprime

    sites = [n == 1 ? ind(gtau[n], 1) : ind(gtau[n], 2) for n in eachindex(gtau)]

    # Apply phase shift to each Qubit
    θ = π * ((-N+1)/N)
    for i in 1:nbit(ft)
        gtau[i] *= noprime(gtau[i] * op("Phase", sites[i]; ϕ= θ * 2^(nbit-i) ))
    end

    # FFT
    return noprime(contract(ft.ftcore.forward, kwargs...))
end