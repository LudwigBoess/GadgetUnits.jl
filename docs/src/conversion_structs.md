```@meta
CurrentModule = GadgetUnits
DocTestSetup = quote
    using GadgetUnits
end
```

# Conversion Structs


GadgetUnits.jl uses Unitful.jl and UnitfulAstro.jl to store the unit conversion factors with actual units in place.
You can convert the internal units of Gadget into cgs units by defining the object [`GadgetPhysicalUnits`](@ref):

```julia
GU = GadgetPhysicalUnits(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;
                         a_scale::T=1.0, hpar::T=1.0,
                         γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T
```

where the keyword arguments are:
- `a_scale::T = 1.0`:  Cosmological scale factor of the simulation. Can be passed with the header `h` as `h.time`.
- `hpar::T = 1.0`:     Hubble constant as 'little h'. Can be passed with header `h` as `h.h0`.
- `γ_th::T = 5.0/3.0`: Adiabatic index of gas.
- `γ_CR::T = 4.0/3.0`: Adiabatic index of cosmic ray component.
- `xH::T = 0.76`:      Hydrogen fraction of the simulation, if run without chemical model.

This returns an object of type [`GadgetPhysicalUnits`](@ref).

If you want to have the same functionality, but without using `Unitful.jl` you can construct a similar object [`GadgetPhysical`](@ref):

```julia
GU = GadgetPhysical(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;
                    a_scale::T=1.0, hpar::T=1.0,
                    γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T
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

```julia
GadgetPhysicalUnits( DT::DataType, 
                     l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;
                     a_scale::Real=1.0, hpar::Real=1.0,
                     γ_th::Real=5.0/3.0, γ_CR::Real=4.0/3.0, xH::Real=0.76)
```

`DT` defaults to `Float64`.

Please be aware that some unit conversions overflow at `Float32`!
