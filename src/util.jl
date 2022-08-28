function siteinds(d::Vector{T}; kwargs...) where {T<:Integer}
    return [siteind(d[n], n; kwargs...) for n in eachindex(d)]
end

function extractsite(x::MPS, n::Int)
    if n == 1
        uniqueind(x[n], x[n+1])
    elseif n == length(x)
        uniqueind(x[n], x[n-1])
    else
        uniqueind(x[n], x[n+1], x[n-1])
    end
end

extractsites(x::MPS) = [extractsite(x, n) for n in eachindex(x)]


function replace_mpo_siteinds!(M::MPO, sites)
    #@show M
    offset(j) = (j == 1) ? 0 : 1
    for j in eachindex(M)
        replaceind!(M[j], ind(M[j], 1+offset(j)), sites[j]')
        replaceind!(M[j], ind(M[j], 2+offset(j)), sites[j])
    end
    return M
end



"""
Reverse the order of the physical indices of a MPO
"""
#revserMPO(reverse([x for x in M]))