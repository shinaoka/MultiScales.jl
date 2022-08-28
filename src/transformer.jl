abstract type AbstractFT end

struct FTCore
    forward::MPO

    function FTCore(nbit)
        sites = siteinds("Qubit", nbit)
        new(_qft(sites))
    end
end

nbit(ft::AbstractFT) = length(ft.ftcore.mpo)