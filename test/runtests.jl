using GadgetUnits, GadgetIO, Test, Unitful, UnitfulAstro, Downloads, Cosmology

Downloads.download("http://www.usm.uni-muenchen.de/~lboess/GadgetIO/snap_144.0", "./snap_144.0")

@testset "GadgetUnits" begin 

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
        end

    end # Conversion Structs

    @testset "Cosmology" begin
        
        h = read_header("snap_144.0")
        h.z = 0.5
        
        @testset "Setup" begin
            c = cosmology(h)
            @test c == Cosmology.FlatLCDM{Float64}(0.704, 0.7279156596296679, 0.272, 8.434037033212373e-5)
        end
        
        @testset "Scale" begin
            c = cosmology(h)

            # tested against script by Klaus Dolag
            @test arcmin_to_kpc(1.0, h)                  ≈ 367.7593389798557
            @test arcmin_to_kpc(c, 1.0, h.z)             ≈ 367.7593389798557
            @test arcmin_to_kpc(1.0u"arcminute", h)      ≈ 367.7593389798557u"kpc"
            @test arcmin_to_kpc(c, 1.0u"arcminute", h.z) ≈ 367.7593389798557u"kpc"
        end

        @testset "Radiation" begin
            c = cosmology(h)

            # tested against script by Klaus Dolag
            @test mJy_to_W(1.0, h)            ≈ 9.681690084092813e23
            @test mJy_to_W(c, 1.0, h.z)       ≈ 9.681690084092813e23
            @test mJy_to_W(1.0u"mJy", h)      ≈ 9.681690084092813e23u"W/Hz"
            @test mJy_to_W(c, 1.0u"mJy", h.z) ≈ 9.681690084092813e23u"W/Hz"
        end

        @testset "Time" begin
            @test age(h)        ≈ 8.69582554048774u"Gyr"
            @test age(h, false) ≈ 8.69582554048774
        end
    end

end

rm("snap_144.0")
