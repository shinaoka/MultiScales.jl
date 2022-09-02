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
    length(linkdims) == length(sites) - 1 || error("Length mismatch")
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