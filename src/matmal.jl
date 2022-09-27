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
           #ab_ = ITensor(ITensors.data(ab_), [sitesa[2n-1], sitesb[2n], newlinks[n]])
           ab_ = ITensor(ITensors.data(ab_), [newsites[n], newlinks[n]])
        elseif n == halfN
           ab_ = permute(ab_, [linksa[2n-2], linksb[2n-2], sitesa[2n-1], sitesb[2n]])
           #ab_ = ITensor(ITensors.data(ab_), [newlinks[n-1], sitesa[2n-1], sitesb[2n]])
           ab_ = ITensor(ITensors.data(ab_), [newlinks[n-1], newsites[n]])
        else
           ab_ = permute(ab_, [linksa[2n-2], linksb[2n-2], sitesa[2n-1], sitesb[2n], linksa[2n], linksb[2n]])
           #ab_ = ITensor(ITensors.data(ab_), [newlinks[n-1], sitesa[2n-1], sitesb[2n], newlinks[n]])
           ab_ = ITensor(ITensors.data(ab_), [newlinks[n-1], newsites[n], newlinks[n]])
        end
        push!(ab_tensors, ab_)
    end
    return splitsiteind(MPS(ab_tensors), siteinds(a))
end