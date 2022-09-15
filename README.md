# Cosmologies

This package calculates basic background cosmology.


## Installation

To install, first press `]` at the REPL to get into the package mode. Then, add the registry and the package,
```julia
pkg> regsitry add https://github.com/Wide-Angle-Team/WATCosmologyJuliaRegistry.git
...
pkg> add Cosmologies
```


## Usage

The modules exposes several kinds of cosmologies. The most basic one is a ΛCDM
cosmology, which is, unsurprisingly, called `ΛCDM`. For example,
```julia
using Cosmologies

cosmo = ΛCDM(h, Ωr, Ωm, Ωk, Ωv)
```
where the arguments have the usual meanings: `h∼0.7` is the Hubble parameter
divided by 100 km/s/Mpc, `Ωr` is the radiation density parameter, `Ωm` the
matter, `Ωk` the curvature, and `Ωv` a cosmological constant.

For convenience, a Planck cosmology is given by
```julia
cosmo = PlanckFlatΛCDM()
```
which will fill in some default values.

Since it is often useful to set the growth rate `f=const` and the growth factor
`D=const`, this module also introduces the `EvolutionlessCosmology`:
```julia
cosmo = EvolutionlessCosmology(cosmobase; D=1, f=1)
```
where `cosmobase` is a regular cosmology that will be used for all functions
except `D` and `f`.

Each cosmology can be used to calculate several quantities. The last character
in the function name indicates the argument, where `a` is the scale factor. The
following functions are exported:
```julia
a = 1 / (1 + z)

Hubble_rate = Hz(cosmo, z)

Omega_radiation_z = Ωra(cosmo, a)
Omega_matter_z = Ωma(cosmo, a)
Omega_curvature_z = Ωka(cosmo, a)
Omega_vaccuum_z = Ωva(cosmo, a)

D = Dz(cosmo, z)
f = fz(cosmo, z)

r = chiz(cosmo, z)  # comoving distance
z = zchi(cosmo, r)  # redshift

comoving_angular_diameter_distance = dAcz(cosmo, z), dLcz
comoving_luminosity_distance = dLcz(cosmo, z)
```
