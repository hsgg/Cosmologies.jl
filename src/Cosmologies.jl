#!/usr/bin/env julia


@doc raw"""
    Cosmologies

This module can be used to calculate basic background cosmology functions for a
ΛCDM cosmology.

Example
=======
```
julia> using Cosmologies
julia> cosmo = PlanckFlatΛCDM()  # alternatively: ΛCDM(h, Ωr, Ωm, Ωk, Ωv)
julia> Dz(cosmo, 0.5)  # linear growth factor D(z)
julia> fz(cosmo, 0.5)  # linear growth rate f(z) = d ln D(z) / d ln a
julia> zchi(cosmo, 1000.0)  # redshift as function of coordinate distance
julia> chiz(cosmo, 0.5)  # coordinate distance
julia> dAcz(cosmo, 0.5)  # comoving angular diameter distance
julia> dAcz(cosmo, 0.5)  # comoving luminosity distance
```
"""
module Cosmologies

export ΛCDM, AbstractΛCDM, AbstractCosmology
export PlanckFlatΛCDM, TimeSliceCosmology, RealSpaceCosmology
export Hz, Ωra, Ωma, Ωka, Ωva, Dz, fz, chiz, zchi, dAcz, dLcz


@deprecate EvolutionlessCosmology(args...; kwargs...) TimeSliceCosmology(args...; kwargs...)


#include("Splines.jl")

#using QuadGK
#using DifferentialEquations
using OrdinaryDiffEq
using Splines

import Base.write
import Base.show


################## constants ##############################
const cm_per_sec = 2.99792458e10  # cm s^-1
const cm_per_Mpc = 3.0857e24      # cm Mpc^-1
const gram_per_Msol = 1.9891e33   # g Msol^-1
const H100 = 100e5 / cm_per_sec   # Mpc^-1

const GNewton = 6.67384e-11 * 1e6 * 1e-3  # cm^3 g^-1 s^-2
const G_Mpc_Msol = GNewton * gram_per_Msol / cm_per_sec^2 / cm_per_Mpc  # Mpc Msol^-1
const a_rad = 7.5657e-15              # erg cm^-3 K^-4 = g cm^-1 s^-2 K^-4
const a_Msol_Mpc3_K4 = a_rad / gram_per_Msol / cm_per_sec^2 * cm_per_Mpc^3


################## cosmologies ############################
abstract type AbstractCosmology end
abstract type AbstractΛCDM <: AbstractCosmology end

# cache object
struct MyCosmologyCache{Tchi,Tglna}
    chi::Tchi
    glna::Tglna
end
function MyCosmologyCache(c::AbstractCosmology)
    return MyCosmologyCache(_make_chi(c), _make_g(c))
end

# for __dot__ syntax:
import Base.length, Base.iterate
length(x::AbstractCosmology) = 1
iterate(x::AbstractCosmology) = x, nothing
iterate(x::AbstractCosmology, i) = nothing

# store cosmological parameters
struct ΛCDM{T<:Real,Tcache,k} <: AbstractΛCDM
    h::T
    Ωr::T
    Ωm::T
    Ωk::T
    Ωv::T
    cache::Tcache
end

ClosedΛCDM{T,Tcache} = ΛCDM{T,Tcache,1}
FlatΛCDM{T,Tcache} = ΛCDM{T,Tcache,0}
OpenΛCDM{T,Tcache} = ΛCDM{T,Tcache,-1}

# default constructors
function ΛCDM(h, Ωr, Ωm, Ωk, Ωv, cache)
    T = eltype(promote(float(h), float(Ωr), float(Ωm), float(Ωk), float(Ωv)))
    Tcache = typeof(cache)
    k = (Ωk == 0) ? 0 : ((Ωk > 0) ? 1 : -1)
    @assert abs(1 - Ωr - Ωm - Ωk - Ωv) < 5eps(T)
    return ΛCDM{T,Tcache,k}(h, Ωr, Ωm, Ωk, Ωv, cache)
end
function ΛCDM(h, Ωr, Ωm, Ωk, Ωv)
    cosmo = ΛCDM(h, Ωr, Ωm, Ωk, Ωv, nothing)
    return ΛCDM(h, Ωr, Ωm, Ωk, Ωv, MyCosmologyCache(cosmo))
end


# Planck cosmology
function PlanckFlatΛCDM(; h=0.6778, Tcmb=2.72548, Ωv=0.69179, Ωnu=3.65e-5)
    #
    ρr0 = a_Msol_Mpc3_K4 * Tcmb^4
    ρc0 = rhocz(ΛCDM(h, 0, 1-Ωv, 0, Ωv, nothing), 0)  # assume Ωr≈0
    Ωr = ρr0 / ρc0 + Ωnu
    Ωm = 1 - Ωv - Ωr
    return ΛCDM(h, Ωr, Ωm, 0, Ωv)
end


############## TimeSliceCosmology
@doc """
    TimeSliceCosmology(basecosmology::AbstractCosmology; D=1, f=1)
    TimeSliceCosmology(h, Ωr, Ωm, Ωk, Ωv; D=1, f=1)

TimeSliceCosmology() is a cosmology where the growth factor `D(z)=const` and
`f(z)=const` for all redshifts. In that sense this describes a static universe.
However, the distance-redshift relation remains without modification from ΛCDM
or whichever is the base cosmology. This is for practical purposes to be able
to define RSD.
"""
struct TimeSliceCosmology{Tcosmo,TD,Tf} <: AbstractCosmology
    cosmobase::Tcosmo
    D::TD
    f::Tf
end

TimeSliceCosmology(; D=1, f=1) = TimeSliceCosmology(PlanckFlatΛCDM(), D, f)

TimeSliceCosmology(cosmobase; D=1, f=1) = TimeSliceCosmology(cosmobase, D, f)

TimeSliceCosmology(h, Ωr, Ωm, Ωk, Ωv; D=1, f=1) = TimeSliceCosmology(ΛCDM(h, Ωr, Ωm, Ωk, Ωv), D, f)


function show(io::IO, c::TimeSliceCosmology)
    print(io, "TimeSliceCosmology(D=$(c.D), f=$(c.f), base=$(c.cosmobase))")
end


# all the functions
for funcname in (:Hz, :rhocz, :Ωra, :Ωma, :Ωka, :Ωva, :chiz, :zchi, :Skz, :zSk, :dAcz, :zdAc, :dLcz)
    @eval $funcname(c::TimeSliceCosmology, x) = $funcname(c.cosmobase, x)
end


# The essential definition of this cosmology:
Dz(c::TimeSliceCosmology, _) = c.D
fz(c::TimeSliceCosmology, _) = c.f


############## RealSpaceCosmology
@doc """
    RealSpaceCosmology(basecosmology::AbstractCosmology; f=1)
    RealSpaceCosmology(h, Ωr, Ωm, Ωk, Ωv; f=1)

`RealSpaceCosmology()` is a cosmology where `f(z)=0`. Everything else is the
same as for a normal cosmology, including the growth factor `D(z)`.
"""
struct RealSpaceCosmology{Tcosmo,Tf} <: AbstractCosmology
    cosmobase::Tcosmo
    f::Tf
end

RealSpaceCosmology(; f=0) = RealSpaceCosmology(PlanckFlatΛCDM(), f)

RealSpaceCosmology(cosmobase; f=0) = RealSpaceCosmology(cosmobase, f)

RealSpaceCosmology(h, Ωr, Ωm, Ωk, Ωv; f=0) = RealSpaceCosmology(ΛCDM(h, Ωr, Ωm, Ωk, Ωv), f)


function show(io::IO, c::RealSpaceCosmology)
    print(io, "RealSpaceCosmology(f=$(c.f), base=$(c.cosmobase))")
end


# all the functions
for funcname in (:Hz, :rhocz, :Ωra, :Ωma, :Ωka, :Ωva, :chiz, :zchi, :Skz, :zSk, :dAcz, :zdAc, :dLcz, :Dz)
    @eval $funcname(c::RealSpaceCosmology, x) = $funcname(c.cosmobase, x)
end


# The essential definition of this cosmology:
fz(c::RealSpaceCosmology, _) = c.f


###################### read/write functions #######################
function write(fname::AbstractString, c::AbstractΛCDM)
    s = ""
    for sym in fieldnames(typeof(c))
        (sym == :cache) && continue
        n = string(sym)
        v = getfield(c, sym)
        s *= "$n = $v\n"
    end
    print(s)
    write(fname, s)
end


function ΛCDM(fname::AbstractString)
    values = Dict()
    open(fname) do f
        for ln in eachline(f)
            s = split(ln, "=")
            length(s) == 1 && continue
            n, v = s
            n = Symbol(strip(n))
            v = strip(v)
            values[n] = try
                v = parse(Float64, v)
            catch
                v
            end
        end
    end
    return ΛCDM(values[:h], values[:Ωr], values[:Ωm], values[:Ωk], values[:Ωv])
end


function show(io::IO, c::ΛCDM)
    print(io, "ΛCDM(h=$(c.h), Ωr=$(c.Ωr), Ωm=$(c.Ωm), Ωk=$(c.Ωk), Ωv=$(c.Ωv))")
end


###################### geometrical functions ######################
# Hubble rate H(z)/h
function Hz(c::AbstractΛCDM, z)
    zp1 = z + 1
    return H100 * sqrt(c.Ωr*zp1^4 + c.Ωm*zp1^3 + c.Ωk*zp1^2 + c.Ωv)
end
function Hz(c::FlatΛCDM, z)
    zp1 = z + 1
    return H100 * sqrt(c.Ωr*zp1^4 + c.Ωm*zp1^3 + c.Ωv)
end


# critical density
function rhocz(c::AbstractΛCDM, z)
    return 3 * c.h^2 * Hz(c, z)^2 / (8π * G_Mpc_Msol)  # Msol Mpc^-3
end


# densities
Ωra(c::AbstractΛCDM, a) = c.Ωr * a^-4 * H100^2 / Hz(c, 1/a-1)^2
Ωma(c::AbstractΛCDM, a) = c.Ωm * a^-3 * H100^2 / Hz(c, 1/a-1)^2
Ωka(c::AbstractΛCDM, a) = c.Ωk * a^-2 * H100^2 / Hz(c, 1/a-1)^2
Ωka(c::FlatΛCDM, a) = 0  # no need for T(0): integers are great!
Ωva(c::AbstractΛCDM, a) = c.Ωv * H100^2 / Hz(c, 1/a-1)^2


# chi
function _make_chi(c::AbstractΛCDM, zmax=1e5)
    # Note: We need an approximate way to handle z<0, because matplotlib needs
    # it when adding redshift as a second x-axis. Here, we simply extrapolate
    # linearly. A proper solution would solve the ODE twice, each time starting
    # at z=0$.
    #I, E = quadgk(z->1/Hz(c,z), 0.0, z)
    f(u,p,t) = 1 / Hz(c, t)
    u0 = 0.0
    tspan = (0.0, zmax)
    prob = ODEProblem(f, u0, tspan)
    sol = solve(prob, Tsit5())
    invsol = Spline1D(sol.u, sol.t, extrapolation=Splines.linear)
    return sol, invsol
end
chiz(c::AbstractΛCDM, z) = c.cache.chi[1](z)
zchi(c::AbstractΛCDM, r) = c.cache.chi[2](r)


# proper distance Sₖ
function Skz(c::ClosedΛCDM, z)
    R0 = 1 / (sqrt(abs(c.Ωk)) * H100)
    return R0 * sin(chiz(c, z) / R0)
end
Skz(c::FlatΛCDM, z) = chiz(c, z)
function Skz(c::OpenΛCDM, z)
    R0 = 1 / (sqrt(abs(c.Ωk)) * H100)
    return R0 * sinh(chiz(c, z) / R0)
end


# inverse of proper distance z(Sₖ)
function zSk(c::ClosedΛCDM, Sk)
    R0 = 1 / (sqrt(abs(c.Ωk)) * H100)
    chi = R0 * asin(Sk / R0)
    return zchi(c, chi)
end
zSk(c::FlatΛCDM, Sk) = zchi(c, Sk)
function zSk(c::OpenΛCDM, Sk)
    R0 = 1 / (sqrt(abs(c.Ωk)) * H100)
    chi = R0 * asinh(Sk / R0)
    return zchi(c, chi)
end


# comoving angular diameter distance
dAcz(c::AbstractΛCDM, z) = Skz(c, z)
zdAc(c::AbstractΛCDM, dAc) = zSk(c, dAc)


# comoving luminosity diameter distance
dLcz(c::AbstractΛCDM, z) = Skz(c, z) * (1 + z)^2


################# growth of structure #####################
function _g_derivs!(du, u, p, t, c::AbstractΛCDM)
        a = exp(t)
        Ωk = Ωka(c, a)
        Ωv = Ωva(c, a)
        weff = -1
        b = 2.5 + 0.5 * (Ωk - 3 * weff * Ωv)
        c = 2 * Ωk + 1.5 * (1 - weff) * Ωv
        g = u[1]
        gp = u[2]
        gpp = -b * gp - c * g
        du[1] = gp
        du[2] = gpp
end
function _make_g(c::AbstractΛCDM)
    f!(du, u, p, t) = _g_derivs!(du, u, p, t, c)

    # from z=30 (MD) to z=0:
    u0 = [1.0, 0.0]
    tspan = (log(1/31), 0.0)  # ln a
    prob = ODEProblem(f!, u0, tspan)
    sol1 = solve(prob, Tsit5(), reltol=0)

    # from z=30 (MD) to z=2000:
    u0 = [1.0, 0.0]
    tspan = (log(1/31), log(1/2001))  # ln a
    prob = ODEProblem(f!, u0, tspan)
    sol2 = solve(prob, Tsit5(), reltol=0)

    u1  = [sol1.u[i][1] for i=1:length(sol1.u)]
    up1 = [sol1.u[i][2] for i=1:length(sol1.u)]
    u2  = [sol2.u[i][1] for i=1:length(sol2.u)]
    up2 = [sol2.u[i][2] for i=1:length(sol2.u)]
    t = [sol1.t; sol2.t]
    idx = unique(i->t[i], sortperm(t))
    t = t[idx]
    u = [u1; u2][idx]
    up = [up1; up2][idx]
    s = Spline1D(t, u)
    sp = Spline1D(t, up)
    return s, sp
end

function Dz(c::AbstractΛCDM, z)
    lna = -log1p(z)
    glna = c.cache.glna
    return exp(lna) * glna[1](lna) / glna[1](0)
end

function fz(c::AbstractΛCDM, z)
    lna = -log1p(z)
    glna = c.cache.glna
    return 1 + glna[2](lna) / glna[1](lna)
end


################# derivatives w.r.t ln(a) #####################
function lnHp(c::FlatΛCDM, lna)
    (-3*c.Ωm)/(2*(c.Ωm + exp(3*lna)*c.Ωv))
end
function lnHpp(c::FlatΛCDM, lna)
    (9*exp(3*lna)*c.Ωm*c.Ωv)/(2*(c.Ωm + exp(3*lna)*c.Ωv)^2)
end
function lnHppp(c::FlatΛCDM, lna)
    (27*exp(3*lna)*c.Ωm*c.Ωv*(c.Ωm - exp(3*lna)*c.Ωv))/(2*(c.Ωm + exp(3*lna)*c.Ωv)^3)
end
function lnHiv(c::FlatΛCDM, lna)
    (81*exp(3*lna)*c.Ωm*c.Ωv*(c.Ωm^2 - 4*exp(3*lna)*c.Ωm*c.Ωv + exp(6*lna)*c.Ωv^2))/
    -  (2*(c.Ωm + exp(3*lna)*c.Ωv)^4)
end


end


# vim: set sw=4 et sts=4 :
