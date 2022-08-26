@doc raw"""
Create a MPO for Fourier transform (in a brute-force way)

The created MPO can transform an input MPS as follows.
We denote the input and output MPS's by ``X`` and ``Y``, respectively.

For inputorder == :normal,

* ``X(x_1, ..., x_N) = X_1(x_1) ... X_N (x_N)``,
* ``Y(y_N, ..., y_1) = Y_1(y_N) ... Y_N (y_1)``.

For inputorder == :reversed,

* ``X(x_N, ..., x_1) = X_1(x_N) ... X_N (x_1)``,
* ``Y(y_1, ..., y_N) = Y_1(y_1) ... Y_N (y_N)``.

Note that the order of the labels of the physical indices is reversed by the transform for smaller bond dimensions of the MPO.

We define two integers using the binary format: ``x = (x_1 x_2 ...., x_N)_2``, ``y = (y_1 y_2 ...., y_N)_2``,
where the right most digits are the least significant digits.

Our definition of the Fourier transform is

```math
    Y(y) = \frac{1}{\sqrt{N}} \sum_{x=0}^{N-1} X(x) e^{s i \frac{2\pi y x}{N}} = \sum_{x=0}^{N-1} T(y, x) X(x),
```

where we define the transformation matrix ``T`` and ``s = \pm 1``.
"""
function _qft(sites; cutoff::Float64=1e-14, sign::Int=1, inputorder=:normal)
    abs(sign) == 1 || error("sign must either 1 or -1")
    inputorder ∈ [:normal, :reversed] || error("Invalid inputourder")

    nbit = length(sites)
    N = 2^nbit

    tmat = zeros(ComplexF64, N, N)
    for t in 0:N-1, x in 0:N-1
        tmat[t+1, x+1] = exp(sign * im * 2π * t * x/N)
    end

    # `tmat`: (y_1, ..., y_N, x_1, ..., x_N)
    tmat ./= sqrt(N)
    tmat = reshape(tmat, ntuple(x->2, 2*nbit))

    if inputorder == :normal
        trans_t = ITensor(tmat, reverse(sites)..., prime(sites)...)
    else
        trans_t = ITensor(tmat, sites..., prime(reverse(sites))...)
    end
    return MPO(trans_t, sites; cutoff=cutoff)
end