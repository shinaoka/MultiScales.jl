function siteinds(d::Vector{T}; kwargs...) where {T<:Integer}
    return [siteind(d[n], n; kwargs...) for n in eachindex(d)]
end