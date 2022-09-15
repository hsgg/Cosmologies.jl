using Cosmologies
using Test

@testset "Cosmologies.jl" begin

    cosmo_spec = [
                  (PlanckFlatΛCDM, ()),
                  (EvolutionlessCosmology, (:D=>0.7, :f=>0.9))
                 ]

    @testset "$cosmology($(join(", ", kwargs)))" for (cosmology, kwargs) in cosmo_spec
        cosmo = cosmology(; kwargs...)

        z = 0.5
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
        @test z ≈ 0.5  atol=1e-4

        comoving_angular_diameter_distance = dAcz(cosmo, z), dLcz
        comoving_luminosity_distance = dLcz(cosmo, z)
    end

end
