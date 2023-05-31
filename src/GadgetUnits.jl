module GadgetUnits

    using Unitful
    using UnitfulAstro
    using Cosmology
    using GadgetIO

    include("conversion_structs/with_units.jl")
    include("conversion_structs/without_units.jl")
    include("cosmo_utility/cosmo.jl")
    include("cosmo_utility/time.jl")
    include("cosmo_utility/radiation.jl")
    include("cosmo_utility/scale.jl")

    export GadgetPhysical,
           GadgetPhysicalUnits,
           cosmology, 
           age, redshift,
           arcmin_to_kpc,
           mJy_to_W


    using SnoopPrecompile    # this is a small dependency

    @precompile_setup begin
        # Putting some things in `setup` can reduce the size of the
        # precompile file and potentially make loading faster.

        h = SnapshotHeader(Int32[13933641, 14181938, 0, 0, 0, 0],
            [0.0, 0.031006304312377084, 0.0, 0.0, 0.0, 0.0],
            0.999000004986117, 0.5, Int32(0), Int32(0),
            UInt32[0xc0000000, 0xc0000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000],
            Int32(0), Int32(2048), 500000.0, 0.307115, 0.692885, 0.6777, Int32(0), Int32(0),
            UInt32[0x00000006, 0x00000006, 0x00000000, 0x00000000, 0x00000000, 0x00000000],
            Int32(0), Int32(0), Int32(3), 0.0f0,
            Int32[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

        @precompile_all_calls begin
            # all calls in this block will be precompiled, regardless of whether
            # they belong to your package or not (on Julia 1.8 and higher)
                
            # conversion structs 
            GadgetPhysical()
            GadgetPhysical(h)
            GadgetPhysicalUnits()
            GadgetPhysicalUnits(h)

            # cosmology 
            c = cosmology(h)
            age(h)
            redshift(12.0, h)
            redshift(12.0u"Gyr", h)

            # scale
            arcmin_to_kpc(0.1, h)
            arcmin_to_kpc(0.1u"°", h)
            arcmin_to_kpc(c, 0.1, 0.5)

            # radiation
            mJy_to_W(100.0, h)
            mJy_to_W(c, 100.0, 0.5)
            mJy_to_W(100.0, 0.5)
        end
    end

end # module
