using Cosmologies
using Test

@testset "Cosmologies.jl" begin

    cosmo_spec = [
                  (PlanckFlatΛCDM, ()),
                  (TimeSliceCosmology, (:D=>0.7, :f=>0.9)),
                  (RealSpaceCosmology, (:f=>0,))
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


@testset "Actual tests" begin

    function test()
        # test construction of PlanckFlatΛCDM()
        c2 = PlanckFlatΛCDM()
        @test typeof(c2) <: Cosmologies.FlatΛCDM

        # read parameters
        h = 0.6778
        Ωr = 1e-5
        Ωk = 0.0
        Ωm = (0.022307 + 0.11865) / h^2
        Ωv = 1 - Ωm - Ωr
        @assert Ωr + Ωm + Ωk + Ωv ≈ 1

        # create cosmology object, simple tests
        c = ΛCDM(h, Ωr, Ωm, Ωk, Ωv)
        @test typeof(c) <: Cosmologies.FlatΛCDM
        @test Hz(c, 0) ≈ Cosmologies.H100
        @test Ωra(c, 1) ≈ Ωr
        @test Ωma(c, 1) ≈ Ωm
        @test Ωka(c, 1) == Ωk == 0
        @test Ωva(c, 1) ≈ Ωv
        @test Dz(c, 0) == 1.0
    end


    function test_io()
        c = PlanckFlatΛCDM()
        mkpath("tmp")
        write("tmp/cosmo.ini", c)

        c2 = ΛCDM("tmp/cosmo.ini")
        @test typeof(c2) <: Cosmologies.FlatΛCDM
        @test Hz(c2, 0) ≈ Cosmologies.H100
        @test Ωka(c2, 1) == 0
        @test Dz(c2, 0) == 1.0
    end

    function test_zchi_chiz()
        h = 0.6777
        Ωr = 0
        Ωm = 0.307115
        Ωk = 0
        Ωv = 1 - Ωr - Ωm - Ωk
        c = ΛCDM(h, Ωr, Ωm, Ωk, Ωv)

        z1 = 0.8
        r1 = chiz(c, z1)
        r2 = 1944.404
        z2 = zchi(c, r2)
        @show z1 z2 r1 r2

        # Are r1,r2 and z1,z2 sorted the same way?
        @test_broken sign(r1-r2) == sign(z1-z2)

        # bigger differences are OK?
        z3 = 0.801
        r3 = chiz(c, z3)
        z4 = 0.799
        r4 = chiz(c, z4)
        @show z3 r3
        @show z4 r4
        @test r3 > r1
        @test r3 > r2
        @test r4 < r1
        @test r4 < r2
    end

    test()
    test_io()
    test_zchi_chiz()

end


@testset "Allocation tests" begin
    c = PlanckFlatΛCDM()
    z = rand(10000)
    @show typeof(c)
    @time chiz(c, z[423])
    @time chiz(c, z[15])
    @time chiz(c, z[1234])
    @time @. chiz(c, z)
    @time @. chiz(c, z)
    @time @. chiz(c, z)
end


# vim: set sw=4 et sts=4 :
