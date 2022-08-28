function siteinds(d::Vector{T}; kwargs...) where {T<:Integer}
    return [siteind(d[n], n; kwargs...) for n in eachindex(d)]
end


extractsites(x::MPS) = [n == 1 ? ind(x[n], 1) : ind(x[n], 2) for n in eachindex(x)]
#extractsites(x::MPS) = [n == 1 ? ind(x[n], 1) : ind(x[n], 2) for n in eachindex(x)]

function replace_mpo_siteinds!(M::MPO, sites)
    #@show M
    offset(j) = (j == 1) ? 0 : 1
    for j in eachindex(M)
        replaceind!(M[j], ind(M[j], 1+offset(j)), sites[j]')
        replaceind!(M[j], ind(M[j], 2+offset(j)), sites[j])
    end
    return M
end