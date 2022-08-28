@doc """
For imaginary-time/-frequency domains
"""
struct ImaginaryTimeFT <: AbstractFT
    ftcore::FTCore

    function ImaginaryTimeFT(ftcore::FTCore)
        new(ftcore)
    end
end


function to_wn(::Fermionic, ft::ImaginaryTimeFT, gtau::MPS, beta::Float64; kwargs...)
    length(gtau) == nbit(ft) || error("Length mismatch")
    nbit_ = length(gtau)
    gtau = noprime(copy(gtau))

    N = 2^nbit_
    sites = extractsites(gtau)

    # Apply phase shift to each Qubit
    θ = π * ((-N+1)/N)
    for i in 1:nbit_
        gtau[i] = noprime(gtau[i] * op("Phase", sites[i]; ϕ= θ * 2^(nbit_-i) ))
    end

    # FFT
    M = forwardmpo(ft.ftcore, sites)
    giv = ITensors.apply(M, gtau; kwargs...)
    giv *= beta * 2^(-nbit_/2)

    return giv
end


function decompose(gtau_smpl::Vector{ComplexF64}, sites; kwargs...)
    nbit = length(sites)
    length(gtau_smpl) == 2^nbit || error("Length mismatch")

    # (g_1, g_2, ...)
    gtau_smpl = reshape(gtau_smpl, repeat([2,], nbit)...)
    gtau_smpl = permutedims(gtau_smpl, reverse(collect(1:nbit)))

    return MPS(gtau_smpl, sites; kwargs...)
end

