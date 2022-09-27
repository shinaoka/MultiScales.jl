function siteinds(d::Vector{T}; kwargs...) where {T<:Integer}
    return [siteind(d[n], n; kwargs...) for n in eachindex(d)]
end

function extractsite(x::Union{MPS,MPO}, n::Int)
    if n == 1
        return noprime(copy(uniqueind(x[n], x[n+1])))
    elseif n == length(x)
        return noprime(copy(uniqueind(x[n], x[n-1])))
    else
        return noprime(copy(uniqueind(x[n], x[n+1], x[n-1])))
    end
end

extractsites(x::Union{MPS,MPO}) = [extractsite(x, n) for n in eachindex(x)]


function replace_mpo_siteinds!(M::MPO, sites_src, sites_dst)
    sites_src = noprime(sites_src)
    sites_dst = noprime(sites_dst)
    for j in eachindex(M)
        replaceind!(M[j], sites_src[j], sites_dst[j])
        replaceind!(M[j], sites_src[j]', sites_dst[j]')
    end
    return M
end



"""
Reverse the order of the physical indices of a MPO
"""
#revserMPO(reverse([x for x in M]))


"""
Create a MPO with ITensor objects of ElType ComplexF64 filled with zero
"""
function _zero_mpo(sites; linkdims=ones(Int, length(sites)-1))
    length(linkdims) == length(sites) - 1 || error("Length mismatch $(length(linkdims)) != $(length(sites)) - 1")
    M = MPO(sites)
    N = length(M)
    links = [Index(1, "n=0,Link")]
    for n in 1:(N-1)
        push!(links, Index(linkdims[n], "n=$(n),Link"))
    end
    push!(links, Index(1, "n=$N,Link"))
    for n in 1:N
        inds_ = (links[n], sites[n]', sites[n], links[n + 1])
        elm_ = zeros(ComplexF64, map(ITensors.dim, inds_)...)
        M[n] = ITensor(elm_, inds_...)
    end
    M[1] *= ITensors.delta(links[1])
    M[N] *= ITensors.delta(links[N + 1])
    
    return M
end


# Compute linkdims for a maximally entangled state
function maxlinkdims(inds)
    N = length(inds)
    for i in 1:N
        @assert !ITensors.hastags(inds, "Link")
    end

    physdims = dim.(inds)

    maxdim = ones(Float64, N-1)
    maxdiml = 1.0
    for i in 1:(N-1)
        maxdiml *= physdims[i]
        maxdim[i] = maxdiml
    end

    maxdimr = 1.0
    for i in 1:(N-1)
        maxdimr *= physdims[N+1-i]
        maxdim[N-i] = min(maxdimr, maxdim[N-i])
    end
    return maxdim
end

linkdims(M) = [dim(ITensors.commonind(M[n], M[n+1])) for n in 1:(length(M)-1)]

links(M) = [ITensors.commonind(M[n], M[n+1]) for n in 1:(length(M)-1)]

function _split(t::ITensor, outerlinks, sites)
    Dleft = dim(outerlinks[1])
    Dright = dim(outerlinks[2])
    prod(size(t)) == 4 * Dleft * Dright || error("Length mismatch")
    sites_ = [outerlinks[1], sites..., outerlinks[2]]
    t_ = ITensor(ITensors.data(t), sites_...)
    U, S, V = svd(t_, sites_[1], sites_[2])
    SV = S * V
    return U, SV
end

function splitsiteind(x::MPS, sites)
    N = length(x)
    2N == length(sites) || error("Length mismatch")
    all((dim(s) for s in siteinds(x)) .== 4) || error("Dim of siteinds must be 4")

    links = [Index(1, "Link=0"), [commonind(x[n], x[n+1]) for n in 1:(N-1)]..., Index(1, "Link=$N")]
    res = ITensor[]
    for n in 1:N
        t1, t2 = _split(x[n], links[n:n+1], sites[2n-1:2n])
        push!(res, t1)
        push!(res, t2)
    end
    linksnew = [links[1], [commonind(res[n], res[n+1]) for n in 1:2N-1]..., links[end]]
    linksnew2 = [Index(dim(linksnew[n+1]), ITensors.defaultlinktags(n)) for n in 0:2N]
    for n in 1:2N
        replaceind!(res[n], linksnew[n], linksnew2[n])
        replaceind!(res[n], linksnew[n+1], linksnew2[n+1])
    end
    res[1] *= onehot(Float64, linksnew2[1]=>1)
    res[end] *= onehot(Float64, linksnew2[end]=>1)
    return MPS(res)
end

function addedges!(x::MPS)
    length(inds(x[1])) == 2 || error("Dim of the first tensor must be 2")
    length(inds(x[end])) == 2 || error("Dim of the last tensor must be 2")
    linkl = Index(1, "Link,l=0")
    linkr = Index(1, "Link,l=$(length(x))")
    x[1] = ITensor(ITensors.data(x[1]), [linkl, inds(x[1])...])
    x[end] = ITensor(ITensors.data(x[end]), [inds(x[end])..., linkr])
    return nothing
end

function removeedges!(x::MPS)
    length(inds(x[1])) == 3 || error("Dim of the first tensor must be 3")
    length(inds(x[end])) == 3 || error("Dim of the last tensor must be 3")
    x[1] *= onehot(Float64, inds(x[1])[1]=>1)
    x[end] *= onehot(Float64, inds(x[end])[end]=>1)
    return nothing
end


function combinesiteinds(x, sites)
    N = length(x)
    halfN = N ÷ 2

    tensors = ITensor[]
    links = [commonind(x[n], x[n+1]) for n in 1:N-1]
    linkdims_new = [dim(links[n]) for n in 2:2:N-1]
    links_new = [
            Index(1, ITensors.defaultlinktags(0)),
            [Index(linkdims_new[n], ITensors.defaultlinktags(n)) for n in 1:(halfN-1)]...,
            Index(1, ITensors.defaultlinktags(halfN))
        ]

    for n in 1:halfN
        t = x[2n-1] * x[2n]
        tdata = Array(t, inds(t)) # Better way?
        push!(tensors, ITensor(tdata, links_new[n], sites[n], links_new[n+1]))
    end

    res =  MPS(tensors)
    removeedges!(res)
    return res
end

_mklinks(dims) = [Index(dims[l], "Link,l=$l") for l in eachindex(dims)]
