| **Documentation**                                                 | **Build Status**                                                                                | **License**                                                                                |
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:| :-----------------------------------------------------------------------------------------------:|
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://LudwigBoess.github.io/GadgetUnits.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LudwigBoess.github.io/GadgetUnits.jl/dev) | [![Build Status](https://travis-ci.org/LudwigBoess/GadgetUnits.jl.svg?branch=master)](https://travis-ci.org/LudwigBoess/GadgetUnits.jl) [![codecov.io](https://codecov.io/gh/LudwigBoess/GadgetUnits.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/GadgetUnits.jl?branch=master) | [![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](LICENSE.md) |

# GadgetUnits.jl

This package is a subproject of [GadJet.jl](https://github.com/LudwigBoess/GadJet.jl) and provides some basic unit conversion functionality to work with the SPH code "Gadget" by Volker Springel (doi:10.1111/j.1365-2966.2005.09655.x).

Unit Conversion
===============

GadgetUnits.jl uses Unitful.jl and UnitfulAstro.jl to store the unit conversion factors with actual units in place.
You can convert the internal units of Gadget into cgs units by defining the object `GadgetPhysicalUnits`:

```julia
GU = GadgetPhysicalUnits(l_unit::Float64=3.085678e21, m_unit::Float64=1.989e43, v_unit::Float64=1.e5;
                         a_scale::Float64=1.0, hpar::Float64=1.0,
                         γ_th::Float64=5.0/3.0, γ_CR::Float64=4.0/3.0, xH::Float64=0.76)
```

where the keyword arguments are:
- `a_scale::Float64 = 1.0`:  Cosmological scale factor of the simulation. Can be passed with the header `h` as `h.time`.
- `hpar::Float64 = 1.0`:     Hubble constant as 'little h'. Can be passed with header `h` as `h.h0`.
- `γ_th::Float64 = 5.0/3.0`: Adiabatic index of gas.
- `γ_CR::Float64 = 4.0/3.0`: Adiabatic index of cosmic ray component.
- `xH::Float64 = 0.76`:      Hydrogen fraction of the simulation, if run without chemical model.

This returns an object of type `GadgetPhysicalUnits` with the following properties:

```julia
struct GadgetPhysicalUnits

    x_cgs::typeof(1.0u"cm")         # position in cm
    v_cgs::typeof(1.0u"cm/s")       # velocity in cm/s
    m_cgs::typeof(1.0u"g")          # mass in g

    t_s::typeof(1.0u"s")            # time in sec
    t_Myr::typeof(1.0u"Myr")        # time in Myr

    E_cgs::typeof(1.0u"erg")        # energy in erg
    E_eV::typeof(1.0u"eV")          # energy in eV

    B_cgs::typeof(1.0u"Gs")         # magnetic field in Gauss

    rho_cgs::typeof(1.0u"g/cm^3")   # density in g/cm^3
    rho_ncm3::typeof(1.0u"n_e")     # density in N_p/cm^3

    T_K::typeof(1.0u"K")            # temperature in K

    P_th_cgs::typeof(1.0u"Ba")      # thermal pressure in Ba
    P_CR_cgs::typeof(1.0u"Ba")      # cosmic ray pressure in Ba

end
```

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


## Primitive unit type

If you want to have the same functionality, but without using `Unitful.jl` you can construct a similar object:

```julia
GU = GadgetPhysical(l_unit::Float64=3.085678e21, m_unit::Float64=1.989e43, v_unit::Float64=1.e5;
                    a_scale::Float64=1.0, hpar::Float64=1.0,
                    γ_th::Float64=5.0/3.0, γ_CR::Float64=4.0/3.0, xH::Float64=0.76)
```

This uses the same conversions, but leaves out the actual unit strings.
