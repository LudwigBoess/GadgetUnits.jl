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
           cosmology, age,
           arcmin_to_kpc,
           mJy_to_W

end # module
