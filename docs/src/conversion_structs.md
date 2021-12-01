```@meta
CurrentModule = GadgetUnits
DocTestSetup = quote
    using GadgetUnits
end
```

# Conversion Structs


GadgetUnits.jl uses Unitful.jl and UnitfulAstro.jl to store the unit conversion factors with actual units in place.
You can convert the internal units of Gadget into cgs units by defining the object [`GadgetPhysicalUnits`](@ref)

```@docs
GadgetPhysicalUnits
```


If you want to have the same functionality, but without using `Unitful.jl` you can construct a similar object [`GadgetPhysical`](@ref):

```@docs
GadgetPhysical
```

This uses the same conversions, but leaves out the actual unit properties.

## Usage

To convert, say positions of gas particles from a cosmological simulation to physical units you can use:

```julia

h     = read_header(filename)

pos   = read_snap(filename, "POS", 0)

GU    = GadgetPhysicalUnits(a_scale=h.time, hpar=h.h0)

pos .*= GU.x_cgs

```

If you have different units than the standard Gadget ones you can call the object cunstructor with different values

```julia
GU = GadgetPhysicalUnits(your_l_unit, your_m_unit, your_v_unit; kwargs...)
```

Converting the units can then be done with Unitful.jl and UnitfulAstro.jl.
So if you want to convert the position units from the default `cm` to `Mpc` you can do this as:

```julia
using Unitful
using UnitfulAstro

pos = read_snap(filename, "POS", 0)
pos = @. pos * GU.x_cgs |> u"Mpc"
```

If you want to get rid of the units, for example if you need basic datatypes again for a function
you can use the funtion `ustrip`:

```julia
pos = ustrip(pos)
```

## Data Types

You can define the datatype of the struct by passing it as an optional first parameter

```@docs
GadgetPhysicalUnits(::DataType)
```

`DT` defaults to `Float64`.

Please be aware that some unit conversions overflow at `Float32`!
