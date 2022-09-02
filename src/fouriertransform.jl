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
    inputorder ∈ [:normal, :reversed] || error("Invalid inputorder")

    nbit = length(sites)
    N = 2^nbit
    sites = noprime(sites)

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
    M = MPO(trans_t, sites; cutoff=cutoff)
    #@show "_qft", M
    return M
end

_delta(i, j) = Int(i == j)

"""
i runs from the left to the right
"""
function _tensor(nbit, layer, i, sign)
    pos_leftmost_phase_tensor = nbit - layer + 1
    if i < pos_leftmost_phase_tensor
        if i == 1
           return reshape(Matrix{ComplexF64}(I, 2, 2), (2, 2, 1))
        else
           return reshape(Matrix{ComplexF64}(I, 2, 2), (1, 2, 2, 1))
        end
    elseif i >= pos_leftmost_phase_tensor
        nphase = i - pos_leftmost_phase_tensor + 1
        ϕ = π * 0.5^(nphase-1)
        _exp(x, k) = exp(sign * im * ϕ * (x-1) * (k-1))

        if i == pos_leftmost_phase_tensor && i == nbit
            arr = zeros(ComplexF64, 1, 2, 2)
            for x in 1:2, k in 1:2
                arr[1,x,k] = _exp(x, k)
            end
            return arr
        end

        if i == pos_leftmost_phase_tensor
            arr = zeros(ComplexF64, 2, 2, 2)
            for x in 1:2, k in 1:2
                arr[x,k,k] = _exp(x, k)
            end
            return arr
        end

        if i == nbit
            # (l, out, in)
            arr = zeros(ComplexF64, 2, 2, 2)
            for x in 1:2, k in 1:2
                arr[k,x,x] = _exp(x, k)
            end
            return arr
        end

        arr = zeros(ComplexF64, 2, 2, 2, 2)
        for x in 1:2, k in 1:2
            arr[k,x,x,k] = _exp(x, k)
        end
        return arr
    end
end



# Note: NOT type stable
function _phasegate(nphase::Int, position, sign=1)
    ϕ = π * 0.5^(nphase-1)
    _exp(x, k) = exp(sign * im * ϕ * (x-1) * (k-1))
    if position == :center
        arr = zeros(ComplexF64, 2, 2, 2, 2)
        for x in 1:2, k in 1:2
            arr[k,x,x,k] = _exp(x, k)
        end
        return arr
    elseif position == :left
        arr = zeros(ComplexF64, 2, 2, 2)
        for x in 1:2, k in 1:2
            arr[x,k,k] = _exp(x, k)
        end
        #@show "left" arr
        return arr
    elseif position == :right
        # (l, out, in)
        arr = zeros(ComplexF64, 2, 2, 2)
        for x in 1:2, k in 1:2
            arr[k,x,x] = _exp(x, k)
        end
        return arr
    elseif position == :only
        arr = zeros(ComplexF64, 1, 2, 2)
        for x in 1:2, k in 1:2
            arr[1,x,k] = _exp(x, k)
        end
        return arr
    else
        error("Invalid position")
    end
end

# Note: NOT type stable
function _identitygate(position)
    arr = Matrix{ComplexF64}(I, 2, 2)
    if position == :center
        return reshape(arr, (1,2,2,1))
    elseif position == :left
        return reshape(arr, (2,2,1))
    else
        error("Invalid position")
    end
end


function _assign!(M::MPO, n::Int, arr; autoreshape=false)
    if autoreshape
        arr = reshape(arr, map(dim, inds(M[n]))...)
    end
    M[n] = ITensor(arr, inds(M[n])...)
    return nothing
end

function _phasegates(sites, ntargetbits, sign)
    N = length(sites)
    offset = length(sites) - ntargetbits
    M = MultiScales._zero_mpo(sites; linkdims=vcat(ones(Int,offset), fill(2,N-1-offset)))

    # I gate
    for n in 1:offset
        pos = (n == 1) ? :left : :center
        _assign!(M, n, _identitygate(pos))
        #if ntargetbits == 1
            #@show "debugI ", offset, M[1]
        #end
    end

    # Phase gates
    if ntargetbits == 1
        _assign!(M, 1 + offset, _phasegate(1, :only, sign))
        #@show "debug ", offset, M[1+offset]
    else
        _assign!(M, 1 + offset, _phasegate(1, :left, sign))
        for i in 2:(ntargetbits-1)
            _assign!(M, i + offset, _phasegate(i, :center, sign))
        end
        _assign!(M, ntargetbits + offset, _phasegate(ntargetbits, :right, sign))
    end
    return M
end

function _qft2(sites; cutoff::Float64=1e-14, sign::Int=1, inputorder=:normal)
    @assert inputorder == :normal
    M = _phasegates(sites, 1, sign)
    for layer in 2:length(sites)
        tmp = _phasegates(sites, layer, sign)
        #M = apply(_phasegates(sites, layer, sign), M; cutoff=cutoff, sign=sign)
        M = apply(tmp, M; cutoff=cutoff, sign=sign)
        for (i, m) in enumerate(tmp)
            println("#############################")
            println("#############################")

            @show i, m
        end
    end
    M *= 2.0^(-0.5 * length(sites))
    return M
end

"""
For length(sites) == 1
The resultant MPO is NOT renormalized.
"""
function _qft_nsite1_wo_norm(sites; sign::Int=1, inputorder=:normal)
    length(sites) == 1 || error("num sites > 1")
    _exp(x, k) = exp(sign * im * π * (x-1) * (k-1))

    arr = zeros(ComplexF64, 2, 2)
    for out in 1:2, in in 1:2
        arr[out, in] = _exp(out, in)
    end

    M = MultiScales._zero_mpo(sites)
    _assign!(M, 1, arr)

    return M
end

function _qft_wo_norm(sites; cutoff::Float64=1e-14, sign::Int=1, inputorder=:normal)
    N = length(sites)
    if N == 1
        return _qft_nsite1_wo_norm(sites; sign=sign, inputorder=inputorder)
    end

    M_prev = _qft_wo_norm(
        sites[2:end]; cutoff=cutoff, sign=sign, inputorder=inputorder)
    M_top = _qft_toplayer(sites; sign=1, inputorder=:normal)

    M = _contract(M_top, M_prev)
    ITensors.truncate!(M; cutoff=cutoff, sign=sign)

    return M
end

function _qft3(sites; cutoff::Float64=1e-14, sign::Int=1, inputorder=:normal)
    M = _qft_wo_norm(sites; cutoff=cutoff, sign=sign, inputorder=inputorder)
    M *= 2.0^(-0.5 * length(sites))
    return M
end


function _qft_toplayer(sites; sign::Int=1, inputorder=:normal)
    N = length(sites)
    N > 1 || error("N must be greater than 1")
    
    tensors = []

    # site = 1
    arr = zeros(ComplexF64, 2, 2, 2)
    for x in 1:2, k in 1:2
        arr[x,k,k] = exp(sign * im * π * (x-1) * (k-1))
    end
    push!(tensors, arr)

    for n in 2:N
        ϕ = π * 0.5^(n-1)
        _exp(x, k) = exp(sign * im * ϕ * (x-1) * (k-1))
        # Right most tensor
        if n == N
            arr = zeros(ComplexF64, 2, 2, 2)
            for x in 1:2, k in 1:2
                arr[k,x,x] = _exp(x, k)
            end
            push!(tensors, arr)
        else
            arr = zeros(ComplexF64, 2, 2, 2, 2)
            for x in 1:2, k in 1:2
                arr[k,x,x,k] = _exp(x, k)
            end
            push!(tensors, arr)
        end
    end

    M = MultiScales._zero_mpo(sites; linkdims=fill(2, N-1))
    for n in 1:N
        _assign!(M, n, tensors[n])
    end

    return M
end



function _contract(M_top, M_prev)
    length(M_top) == length(M_prev) + 1 || error("Length mismatch")
    M_top = ITensors.replaceprime(M_top, 1=>2; tags="Qubit")
    M_top = ITensors.replaceprime(M_top, 0=>1; tags="Qubit")
    M_top_ = ITensors.data(M_top)
    M_prev_ = ITensors.data(M_prev)

    M_data = [M_top_[1]]
    for n in 1:length(M_prev_)
        push!(M_data, M_top_[n+1] * M_prev_[n])
    end

    M =  MPO(M_data)
    M =  ITensors.replaceprime(M, 1=>0; tags="Qubit")
    M =  ITensors.replaceprime(M, 2=>1; tags="Qubit")

    return M
end