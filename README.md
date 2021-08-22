| **Documentation**                                                 | **Build Status**                                                                                | **License**                                                                                |
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:| :-----------------------------------------------------------------------------------------------:|
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://LudwigBoess.github.io/GadgetUnits.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LudwigBoess.github.io/GadgetUnits.jl/dev) | [![Run CI on master](https://github.com/LudwigBoess/GadgetUnits.jl/actions/workflows/jlpkgbutler-ci-master-workflow.yml/badge.svg)](https://github.com/LudwigBoess/GadgetUnits.jl/actions/workflows/jlpkgbutler-ci-master-workflow.yml) [![codecov.io](https://codecov.io/gh/LudwigBoess/GadgetUnits.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/GadgetUnits.jl?branch=master) | [![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](LICENSE.md) |

# GadgetUnits.jl

This package is a subproject of [GadJet.jl](https://github.com/LudwigBoess/GadJet.jl) and provides some basic unit conversion functionality to work with the SPH code "Gadget" by Volker Springel (doi:10.1111/j.1365-2966.2005.09655.x).

Unit Conversion
===============

GadgetUnits.jl uses Unitful.jl and UnitfulAstro.jl to store the unit conversion factors with actual units in place.
You can convert the internal units of Gadget into cgs units by defining the object `GadgetPhysicalUnits`:

```julia
GU = GadgetPhysicalUnits(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;
                         a_scale::T=1.0, hpar::T=1.0,
                         γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T
```

where the keyword arguments are:
- `a_scale = 1.0`:  Cosmological scale factor of the simulation. Can be passed with the header `h` as `h.time`.
- `hpar = 1.0`:     Hubble constant as 'little h'. Can be passed with header `h` as `h.h0`.
- `γ_th = 5.0/3.0`: Adiabatic index of gas.
- `γ_CR = 4.0/3.0`: Adiabatic index of cosmic ray component.
- `xH = 0.76`:      Hydrogen fraction of the simulation, if run without chemical model.

This returns an object of type `GadgetPhysicalUnits` with the following properties:

| Name                     | Meaning                        |
|: ----------------------- |:-----------------------------  |
| `x_cgs::T`         | position in cm                 |
| `x_kpc::T`         | position in kpc                |
| `v_cgs::T`         | velocity in cm/s               |
| `v_kms::T`         | velocity in km/s               |
| `m_cgs::T`         | mass in g                      |
| `m_msun::T`        | mass in Msun                   |
| `t_s::T`           | time in sec                    |
| `t_Myr::T`         | time in Myr                    |
| `E_cgs::T`         | energy in erg                  |
| `E_eV::T`          | energy in eV                   |
| `B_cgs::T`         | magnetic field in Gauss        |
| `rho_cgs::T`       | density in ``g/cm^3``          |
| `rho_ncm3::T`      | density in ``n_p/cm^3``        |
| `T_K::T`           | temperature in K               |
| `T_eV::T`          | temperature in eV              |
| `P_th_cgs::T`      | thermal pressure in Ba         |
| `P_CR_cgs::T`      | cosmic ray pressure in Ba      |

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
GU = GadgetPhysical(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;
                    a_scale::T=1.0, hpar::T=1.0,
                    γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T
```

This uses the same conversions, but leaves out the actual unit strings.
