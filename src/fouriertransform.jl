@doc raw"""
Create a MPO for QFT in a brute force way

F(t) = \sqrt{N} \sum_{x=0}^{N-1} f(x) e^{s * i \frac{2\pi t x}{N}} = \sum_{x=0}^{N-1} T(t, x) f(x),

where N = 2^nbit and `sign` defaults to 1.

In the input `x`, the first bit is to the most significant digit.
In the output `t`, the order of the qubits is reversed.
"""
function _qft(sites; cutoff::Float64=1e-14, sign::Int=1)
    nbit = length(sites)
    N = 2^nbit

    tmat = zeros(ComplexF64, N, N)
    for t in 0:N-1, x in 0:N-1
        tmat[t+1, x+1] = exp(sign * im * 2π * t * x/N)
    end

    # `tmat`: (t_0, ..., t_{Q-1}, x_0, ..., x_{Q-1})
    tmat ./= sqrt(N)
    tmat = reshape(tmat, ntuple(x->2, 2*nbit))

    trans_t = ITensor(tmat, reverse(sites)..., prime(sites)...)
    return MPO(trans_t, sites; cutoff=cutoff)
end