using GadgetUnits, GadgetIO, Test, Unitful, UnitfulAstro, Cosmology


@testset "GadgetUnits" begin 

    # define header here
    h = SnapshotHeader(Int32[13933641, 14181938, 0, 0, 0, 0],
        [0.0, 0.031006304312377084, 0.0, 0.0, 0.0, 0.0],
        0.999000004986117, 0.5, Int32(0), Int32(0),
        UInt32[0xc0000000, 0xc0000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000],
        Int32(0), Int32(2048), 500000.0, 0.307115, 0.692885, 0.6777, Int32(0), Int32(0),
        UInt32[0x00000006, 0x00000006, 0x00000000, 0x00000000, 0x00000000, 0x00000000],
        Int32(0), Int32(0), Int32(3), 0.0f0,
        Int32[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    @testset "Unit Structs" begin

        @testset "GadgetPhysicalUnits" begin

            @testset "Float64" begin 
                GU = GadgetPhysicalUnits(Float64)

                @test GU.t_s ≈ 3.085678e16u"s"
                @test GU.E_cgs ≈ 1.989e53u"erg"
            end

            @testset "Float32" begin 
                GU = GadgetPhysicalUnits(Float32)
                @test GU.t_s        ≈ 3.0856778f16u"s"
            end

            @testset "SnapshotHeader" begin

                @testset "cosmo" begin
                    h.z  = 0.5
                    GU = GadgetPhysicalUnits(h)
                    @test GU.x_cgs ≈ 3.0354414441001423e21u"cm"
                end

                @testset "no cosmo" begin
                    h.z  = 0.0
                    GU = GadgetPhysicalUnits(h)
                    @test GU.t_s ≈ 4.553162166150214e16u"s"
                end
                
            end
        end

        
        @testset "GadgetPhysical" begin

            @testset "Float64" begin 
                GU = GadgetPhysical(Float64)

                @test GU.t_s   ≈ 3.085678e16
                @test GU.E_cgs ≈ 1.989e53
            end

            @testset "Float32" begin 
                GU = GadgetPhysical(Float32)

                @test GU.t_s        ≈ 3.085678e16
                @test GU.x_physical ≈ 1.0
            end

            @testset "SnapshotHeader" begin

                @testset "cosmo" begin
                    h.z = 0.5
                    GU = GadgetPhysical(h)
                    @test GU.x_cgs ≈ 3.0354414441001423e21
                end

                @testset "no cosmo" begin
                    h.z = 0.0
                    GU = GadgetPhysical(h)
                    @test GU.t_s        ≈ 4.553162166150214e16
                end
                
            end
        end

    end # Conversion Structs

    @testset "Cosmology" begin
        
        h.z = 0.5
        
        @testset "Setup" begin
            c = cosmology(h)
            @test c.h     ≈ h.h0
            @test c.Ω_m   ≈ h.omega_0
            #Cosmology.FlatLCDM{Float64}(0.704, 0.7279156596296679, 0.272, 8.434037033212373e-5)
        end
        
        @testset "Scale" begin
            c = cosmology(h)

            # tested against script by Klaus Dolag
            @test arcmin_to_kpc(1.0, h)                  ≈ 377.3573133756037
            @test arcmin_to_kpc(c, 1.0, h.z)             ≈ 377.3573133756037
            @test arcmin_to_kpc(1.0u"arcminute", h)      ≈ 377.3573133756037u"kpc"
            @test arcmin_to_kpc(c, 1.0u"arcminute", h.z) ≈ 377.3573133756037u"kpc"
        end

        @testset "Radiation" begin
            c = cosmology(h)
            
            # Luminsity Distance in Mpc
            dL = luminosity_dist(c, h.z)

            # reference values tested against script by Klaus Dolag
            @test mJy_to_W(1.0, h)       ≈ 1.0193640161171312e24
            @test mJy_to_W(1.0u"mJy", h) ≈ 1.0193640161171312e24u"W/Hz"

            # cross-check
            # without untis
            @test mJy_to_W(c, 1.0, h.z)           ≈ mJy_to_W(1.0, h)
            @test mJy_to_W(1.0, ustrip(dL))       ≈ mJy_to_W(1.0, h)
            @test mJy_to_W(1.0u"mJy", ustrip(dL)) ≈ mJy_to_W(1.0, h)

            # with units
            @test mJy_to_W(c, 1.0u"mJy", h.z) ≈ mJy_to_W(1.0u"mJy", h)
            @test mJy_to_W(1.0u"mJy", dL)     ≈ mJy_to_W(1.0u"mJy", h)
        end

        @testset "Time" begin
            @test age(h)        ≈ 8.61908274173147u"Gyr"
            @test age(h, false) ≈ 8.61908274173147

            @test redshift(8.61908274173147u"Gyr", h) ≈ 0.5
            @test redshift(8.61908274173147, h)       ≈ 0.5

            c = cosmology(h)
            @test redshift(8.61908274173147u"Gyr", c) ≈ 0.5
            @test redshift(8.61908274173147, c)       ≈ 0.5
        end
    end

end
