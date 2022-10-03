function matmul(a::MPS, b::MPS)
    N = length(a)
    mod(N, 2) == 0 || error("Length of a must be even")
    length(a) == length(b) || error("Length mismatch")
    halfN = N ÷ 2

    ab_tensors = ITensor[]
    sitesa = siteinds(a)
    sitesb = siteinds(b)
    linksa = MultiScales.links(a)
    linksb = MultiScales.links(b)
    newlinks = _mklinks([dim(linksa[2n]) * dim(linksb[2n]) for n in 1:halfN-1])
    newsites = [Index(4, "csite=$s") for s in 1:halfN]
    
    for n in 1:halfN
        a_ = a[2*n-1:2*n] # copy
        b_ = b[2*n-1:2*n]
        cind = Index(2, "Qubit")
        replaceind!(a_[2], sitesa[2*n], cind)
        replaceind!(b_[1], sitesb[2*n-1], cind)
        ab_ = a_[1] * ((a_[2] * b_[1]) * b_[2])
        if n == 1
           ab_ = permute(ab_, [sitesa[2n-1], sitesb[2n], linksa[2n], linksb[2n]])
           ab_ = ITensor(ITensors.data(ab_), [newsites[n], newlinks[n]])
        elseif n == halfN
           ab_ = permute(ab_, [linksa[2n-2], linksb[2n-2], sitesa[2n-1], sitesb[2n]])
           ab_ = ITensor(ITensors.data(ab_), [newlinks[n-1], newsites[n]])
        else
           ab_ = permute(ab_, [linksa[2n-2], linksb[2n-2], sitesa[2n-1], sitesb[2n], linksa[2n], linksb[2n]])
           ab_ = ITensor(ITensors.data(ab_), [newlinks[n-1], newsites[n], newlinks[n]])
        end
        push!(ab_tensors, ab_)
    end
    return splitsiteind(MPS(ab_tensors), siteinds(a))
end

function _tompo_matmul(t1::ITensor, t2::ITensor, sites, links, s)
    s1, s2 = sites
    l0, l1, l2 = links

    s1_new = Index(2, "site1")
    s2_new = Index(2, "site2")

    t1_ = copy(t1)
    replaceind!(t1_, s1, s1_new')

    t2_ = copy(t2)
    replaceind!(t2_, s2, s1_new)

    res = (t1_ * t2_) * delta(s2_new', s2_new)

    res = permute(res, [l0, l2, s1_new', s2_new', s1_new, s2_new])
    res = ITensor(ITensors.data(res), [l0, l2, s', s])
    return res
end


function tompo_matmul(a::MPS, csites)
    N = length(a)
    sites = siteinds(a)

    mod(N, 2) == 0 || error("Length of a must be even")
    length(a) == 2*length(csites) || error("Length mismatch")
    halfN = N ÷ 2

    addedges!(a)
    linksa = _linkinds(a, sites)
    tensors = [
                _tompo_matmul(
                    a[2*n-1], a[2*n], sites[2*n-1:2*n], linksa[2*n-1:2*n+1], csites[n])
                for n in 1:halfN
            ]
    #for (n, t) in enumerate(tensors)
        #@show n, t
    #end
    removeedges!(a, sites)

    res = MPO(tensors)
    removeedges!(res, csites)
    return res
end