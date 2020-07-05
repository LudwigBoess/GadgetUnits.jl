using GadgetUnits, Test, Unitful

@testset "Unit Conversion" begin

    GU = GadgetPhysicalUnits()

    @test GU.t_s ≈ 3.085678e16u"s"
    @test GU.E_cgs ≈ 1.989e53u"erg"

    d = strip_unit(1.0u"g")
    @test d == 1.0

    GU = GadgetPhysical()

    @test GU.t_s ≈ 3.085678e16
    @test GU.E_cgs ≈ 1.989e53


end