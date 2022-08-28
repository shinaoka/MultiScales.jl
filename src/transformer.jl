abstract type AbstractFT end

struct FTCore
    forward::MPO

    function FTCore(sites)
        new(_qft(sites))
    end
end

nbit(ft::AbstractFT) = length(ft.ftcore.forward)


function forwardmpo(ftcore::FTCore, sites)
    M = copy(ftcore.forward)
    replace_mpo_siteinds!(M, sites)
    return M
end