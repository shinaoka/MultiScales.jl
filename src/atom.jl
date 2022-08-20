module atom
import ..SparseBSE: check_fermionic, check_ph_convention, asarray
using OMEinsum: OMEinsum
import DoubleFloats: Double64

_delta(n::Int64, npr::Int64)::Int64 = (n == npr) ? 1 : 0

# TODO
# * Flip the sign of wb
# * Exchage v and v + w in X^ph_0

"""
Return a single line in the atomic limit
"""
function gf(U::Float64, beta::Float64, n::Int64)
    n = check_fermionic(n)
    Uhalf = 0.5 * U
    nu = π / beta * n
    return -im * nu / (nu^2 + Uhalf^2)
end

"""
Return the self-energy
"""
function sigma(U::Float64, beta::Float64, n::Int64)::ComplexF64
    n = check_fermionic(n)
    Uhalf = 0.5 * U 
    nu = π / beta * n 
    return -im*(Uhalf^2)/nu + Uhalf
end

function sigma_spin_basis(U::Float64, beta::Float64, n::Int64)::Matrix{ComplexF64}
    res = zeros(ComplexF64, 2, 2)
    res[1,1] = res[2,2] = sigma(U, beta, n)
    return res
end

function gf_spin_basis(U::Float64, beta::Float64, n::Int64)::Matrix{ComplexF64}
    res = zeros(ComplexF64, 2, 2)
    res[1, 1] = res[2, 2] = gf(U, beta, n)
    return res
end

"""
Irreducible vertex in the density/magnetic repr and the ph convention.

Implements Eq. 19 of Phys. Rev. B 98, 235107 (2018).
"""
function Γ_ph(channel::Symbol, U_::Float64, beta_::Float64,
              n::Int64, npr::Int64, m::Int64)::ComplexF64
    U = Double64(U_)
    beta = Double64(beta_)
    n, npr, m = check_ph_convention(n, npr, m)
    nu = n * (π / beta)
    nupr = npr * (π / beta)
    w = m * (π / beta)

    uhalf2 = (U / 2)^2
    expbu = exp(beta * U / 2)

    if channel === :d
        ar2 = 3 * uhalf2
        br2 = uhalf2 * (-1 + 3 * expbu) / (1 + expbu)
        a0r = 1
        b0r = 1
    elseif channel === :m
        ar2 = -uhalf2
        br2 = uhalf2 * (-expbu + 3) / (1 + expbu)
        a0r = 1
        b0r = 1
    else
        error("Invalid channel" * channel)
    end

    sign::Double64 = Dict(:d => 1, :m => -1)[channel]

    result::Complex{Double64} = sign * U

    sqrt_w2 = sqrt(Complex{Double64}(4 * br2 + w^2))
    @assert !isnan(sqrt_w2)
    result -= U * uhalf2 * (uhalf2 * (br2 / uhalf2 + 1)^2 + w^2) /
              ((U * tan(beta / 4 * (sqrt_w2 + w))) / sqrt_w2 + sign) /
              (nu * (nu + w) - br2) /
              (nupr * (nupr + w) - br2)
    @assert !isnan(result)

    enum = (beta / 2) * (nu^2 + uhalf2) * ((nu + w)^2 + uhalf2)
    term_a = ar2 * enum / (a0r * (nu * (nu + w) - ar2) * nu * (nu + w))
    term_b = br2 * enum / (b0r * (nu * (nu + w) - br2) * nu * (nu + w))
    @assert !isnan(term_a)
    @assert !isnan(term_b)

    result += (n == npr ? term_b + term_a : 0.0)
    result += (n == -npr - m ? term_b - term_a : 0.0)

    return result
end

function Γ_ph_uu_ud(U::Float64, beta::Float64,
                    n::Int64, npr::Int64, m::Int64)::Tuple{ComplexF64,ComplexF64}
    d = Γ_ph(:d, U, beta, n, npr, m)
    m = Γ_ph(:m, U, beta, n, npr, m)
    return 0.5 * (d + m), 0.5 * (d - m)
end

"""
Full vertex in the density/magnetic repr and the ph convention.

Implements Eq. 27 of Phys. Rev. B 98, 235107 (2018).
"""
function F_ph_uu_ud(U, beta, n, npr, m)::Tuple{ComplexF64,ComplexF64}
    n, npr, m = check_ph_convention(n, npr, m)
    nu = n * (π / beta)
    nupr = npr * (π / beta)
    w = m * (π / beta)

    uhalf2 = (U / 2)^2
    expbu = exp(beta * U / 2)
    denom = nu * (nu + w) * nupr * (nupr + w)

    # up-up component
    Fuu = beta * uhalf2 * (nu^2 + uhalf2) * ((nupr + w)^2 + uhalf2) /
          (nu^2 * (nupr + w)^2) * (_delta(n, npr) - _delta(m, 0))

    # up-down component
    Fud = U
    Fud -= 1 / 8 * U^3 * (nu^2 + (nu + w)^2 + nupr^2 + (nupr + w)^2) / denom
    Fud -= 3 / 16 * U^5 / denom
    Fud -= beta * uhalf2 / (1 + expbu) *
           ((nu + w)^2 + uhalf2) * ((nupr + w)^2 + uhalf2) /
           ((nu + w)^2 * (nupr + w)^2) * (2 * _delta(n, -npr - m) + _delta(m, 0))
    Fud += beta * uhalf2 / (1 + 1 / expbu) *
           (nu^2 + uhalf2) * ((nupr + w)^2 + uhalf2) /
           (nu^2 * (nupr + w)^2) * (2 * _delta(n, npr) + _delta(m, 0))
    return Fuu, Fud
end

function F_ph(channel::Symbol, U, beta, n, npr, m)::ComplexF64
    Fuu, Fud = F_ph_uu_ud(U, beta, n, npr, m)
    return Dict(:d => Fuu + Fud, :m => Fuu - Fud)[channel]
end

# Definition of functions F_ph_spin_basis, Γ_ph_spin_basis
function F_ph_spin_basis(U::Float64, beta::Float64, n::Int64, npr::Int64,
                         m::Int64)::Array{ComplexF64,4}
    uu, ud = F_ph_uu_ud(U, beta, n, npr, m)
    barud = F_ph_uu_ud(U, beta, n, n + m, npr - n)[2]
    res = zeros(ComplexF64, 2, 2, 2, 2)
    res[1, 1, 1, 1] = res[2, 2, 2, 2] = uu
    res[1, 1, 2, 2] = res[2, 2, 1, 1] = ud
    res[2, 1, 1, 2] = res[1, 2, 2, 1] = -barud
    return res
end

function Γ_ph_spin_basis(U::Float64, beta::Float64, n::Int64, npr::Int64,
                         m::Int64)::Array{ComplexF64,4}
    uu, ud = Γ_ph_uu_ud(U, beta, n, npr, m)
    res = zeros(ComplexF64, 2, 2, 2, 2)
    res[1, 1, 1, 1] = res[2, 2, 2, 2] = uu
    res[1, 1, 2, 2] = res[2, 2, 1, 1] = ud
    res[2, 1, 1, 2] = res[1, 2, 2, 1] = uu - ud
    return res
end

""" Bare suspectibility """
function X0_ph(U::Float64, beta::Float64, n::Int64, npr::Int64, m::Int64)::ComplexF64
    n, npr, m = check_ph_convention(n, npr, m)
    gfat(n) = gf(U, beta, n)
    return n == npr ? beta * gfat(n) * gfat(n + m) : 0.0
end

function X0_ph_spin_basis(U::Float64, beta::Float64, n::Int64, npr::Int64,
                          m::Int64)::Array{ComplexF64,4}
    res = zeros(ComplexF64, 2, 2, 2, 2)
    for i in 1:2, j in 1:2
        res[i, j, j, i] = X0_ph(U, beta, n, npr, m)
    end
    return res
end

""" Disconnected part of Four-point Green's function """
function g4pt_disconn_ph_spin_basis(U::Float64, beta::Float64, n::Int64, npr::Int64,
                                    m::Int64)::Array{ComplexF64,4}
    g1 = gf_spin_basis(U, beta, n + m)
    g2 = gf_spin_basis(U, beta, npr)
    OMEinsum.@ein term1[a, b, c, d] := g1[a, b] * g2[c, d]
    OMEinsum.@ein term2[a, b, c, d] := g1[a, d] * g2[c, b]
    return beta * (term1 * _delta(m, 0) - term2 * _delta(n, npr))
end

""" Four-point Green's function """
function g4pt_ph_spin_basis(U::Float64, beta::Float64, n::Int64, npr::Int64,
                            m::Int64)::Array{ComplexF64,4}
    term1 = g4pt_disconn_ph_spin_basis(U, beta, n, npr, m)

    X_L = X0_ph_spin_basis(U, beta, n, n, m)
    X_R = X0_ph_spin_basis(U, beta, npr, npr, m)
    F = F_ph_spin_basis(U, beta, n, npr, m)
    OMEinsum.@ein term2[a, b, c, d] := X_L[b, a, e, f] * F[e, f, g, h] * X_R[g, h, d, c]

    return term1 - (beta^(-2)) * term2
end

"""
Generalized suspectibility
    computed from full vertex F

X = X0 + X^conn,
    where X^conn = (1/β^2) * X0 * F * X0.
"""
function X_ph_spin_basis(U::Float64, beta::Float64, n::Int64, npr::Int64,
                         m::Int64)::Array{ComplexF64,4}
    X0 = X0_ph_spin_basis(U, beta, n, npr, m)

    X0_L = X0_ph_spin_basis(U, beta, n, n, m)
    F = F_ph_spin_basis(U, beta, n, npr, m)
    X0_R = X0_ph_spin_basis(U, beta, npr, npr, m)
    OMEinsum.@ein X_conn[a, b, c, d] := X0_L[a, b, e, f] * F[e, f, g, h] * X0_R[g, h, c, d]
    X_conn /= beta^2

    return X0 + X_conn
end

# Vectorized version
for func in [:Γ_ph_spin_basis, :X_ph_spin_basis, :F_ph_spin_basis, :X0_ph_spin_basis,
             :g4pt_ph_spin_basis, :g4pt_disconn_ph_spin_basis]
    @eval begin
        function $(func)(U::Float64, beta::Float64, n::Vector{Int64}, npr::Vector{Int64},
                         m::Vector{Int64})::Array{ComplexF64,5}
            return asarray($(func).(U, beta, n, npr, m))
        end
    end
end

end
