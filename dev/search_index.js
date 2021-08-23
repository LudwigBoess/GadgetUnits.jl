var documenterSearchIndex = {"docs":
[{"location":"api/#API-Reference","page":"API reference","title":"API Reference","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"CurrentModule = GadgetUnits\nDocTestSetup = quote\n    using GadgetUnits\nend","category":"page"},{"location":"api/","page":"API reference","title":"API reference","text":"","category":"page"},{"location":"api/#Exported-Functions","page":"API reference","title":"Exported Functions","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [GadgetUnits]\nPrivate = false\nOrder = [:function]","category":"page"},{"location":"api/#Cosmology.age","page":"API reference","title":"Cosmology.age","text":"age(h::SnapshotHeader)\n\nComputes the age of the universe given the properties from SnapshotHeader.\n\n\n\n\n\n","category":"function"},{"location":"api/#Cosmology.cosmology-Tuple{GadgetIO.SnapshotHeader}","page":"API reference","title":"Cosmology.cosmology","text":"cosmology(h::SnapshotHeader)\n\nDefines a Cosmology.AbstractCosmology from the properties of h.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.arcmin_to_kpc-Tuple{Cosmology.AbstractCosmology, Real, Real}","page":"API reference","title":"GadgetUnits.arcmin_to_kpc","text":"arcmin_to_kpc(c::Cosmology.AbstractCosmology, θ::Real, z::Real)\n\nConvert arcmin to kpc for a given redshift z and cosmology c. Deletes unit information.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.arcmin_to_kpc-Tuple{Cosmology.AbstractCosmology, Unitful.AbstractQuantity, Real}","page":"API reference","title":"GadgetUnits.arcmin_to_kpc","text":"arcmin_to_kpc(c::Cosmology.AbstractCosmology, θ::Unitful.AbstractQuantity, z::Real)\n\nConvert arcmin to kpc for a given redshift z and cosmology c. Conserves unit information.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.arcmin_to_kpc-Tuple{Real, GadgetIO.SnapshotHeader}","page":"API reference","title":"GadgetUnits.arcmin_to_kpc","text":"arcmin_to_kpc(θ::Real, h::SnapshotHeader)\n\nConvert arcmin to kpc for a given redshift and cosmology taken from SnapshotHeader.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.arcmin_to_kpc-Tuple{Unitful.AbstractQuantity, GadgetIO.SnapshotHeader}","page":"API reference","title":"GadgetUnits.arcmin_to_kpc","text":"arcmin_to_kpc(θ::Unitful.AbstractQuantity, h::SnapshotHeader)\n\nConvert arcmin to kpc for a given redshift and cosmology taken from SnapshotHeader.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.mJy_to_W-Tuple{Cosmology.AbstractCosmology, Real, Real}","page":"API reference","title":"GadgetUnits.mJy_to_W","text":"mJy_to_W( c::Cosmology.AbstractCosmology, S::Real, z::Real )\n\nConverts synchrotron emission S in [mJy] to [W/Hz] for a given redshift z and cosmology c. Deletes unit information.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.mJy_to_W-Tuple{Cosmology.AbstractCosmology, Unitful.AbstractQuantity, Real}","page":"API reference","title":"GadgetUnits.mJy_to_W","text":"mJy_to_W( c::Cosmology.AbstractCosmology, S::Unitful.AbstractQuantity, z::Real )\n\nConverts synchrotron emission S in [mJy] to [W/Hz] for a given redshift z and cosmology c. Conserves unit information.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.mJy_to_W-Tuple{Union{Real, Unitful.AbstractQuantity}, GadgetIO.SnapshotHeader}","page":"API reference","title":"GadgetUnits.mJy_to_W","text":"mJy_to_W(S::Union{Real, Unitful.AbstractQuantity}, h::SnapshotHeader)\n\nConverts synchrotron emission S in [mJy] to [W/Hz] for a given redshift and cosmology taken from SnapshotHeader.\n\n\n\n\n\n","category":"method"},{"location":"api/#Exported-Types","page":"API reference","title":"Exported Types","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [GadgetUnits]\nPrivate = false\nOrder = [:type]","category":"page"},{"location":"api/#GadgetUnits.GadgetPhysical","page":"API reference","title":"GadgetUnits.GadgetPhysical","text":"GadgetPhysical( DT::DataType, \n                l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;\n                a_scale::Real=1.0, hpar::Real=1.0,\n                γ_th::Real=5.0/3.0, γ_CR::Real=4.0/3.0, xH::Real=0.76)\n\nSet up a unit struct with a given DataType.\n\n\n\n\n\n","category":"type"},{"location":"api/#GadgetUnits.GadgetPhysical-Union{Tuple{}, Tuple{T}, Tuple{T, T}, Tuple{T, T, T}} where T","page":"API reference","title":"GadgetUnits.GadgetPhysical","text":"GadgetPhysical(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;\n            a_scale::T=1.0, hpar::T=1.0,\n            γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T\n\nCreates a datatype GadgetPhysical which holds the conversion factors between comoving code units and physical units, without unit information.\n\nExamples\n\njulia> GU = GadgetPhysical()\nGadgetPhysicalUnits(3.085678e21, 100000.0, 1.989e43, 3.085678e16, 977.7923542981722, 1.989e53, 1.2414361549102458e65, 1.0, 6.769911178294544e-22, 743179.9340255889, 47.50882854026919, 6.769911178294544e-12, 6.769911178294544e-12)\njulia> GU.x_cgs\n3.085678e21\n\nKeyword arugments specify:\n\nArguments\n\na_scale::T = 1.0:  Cosmological scale factor of the simulation. Can be passed with the header h as h.time.\nhpar::T = 1.0:     Hubble constant as 'little h'. Can be passed with header h as h.h0.\nγ_th::T = 5.0/3.0: Adiabatic index of gas.\nγ_CR::T = 4.0/3.0: Adiabatic index of cosmic ray component.\nxH::T = 0.76:      Hydrogen fraction of the simulation, if run without chemical model.\n\nFields\n\n| Name                     | Meaning                        | |: –––––––––––- |:––––––––––––––-  | | x_cgs::T         | position in cm                 | | x_kpc::T         | position in kpc                | | v_cgs::T         | velocity in cm/s               | | v_kms::T         | velocity in km/s               | | m_cgs::T         | mass in g                      | | m_msun::T        | mass in Msun                   | | t_s::T           | time in sec                    | | t_Myr::T         | time in Myr                    | | E_cgs::T         | energy in erg                  | | E_eV::T          | energy in eV                   | | B_cgs::T         | magnetic field in Gauss        | | rho_cgs::T       | density in gcm^3          | | rho_ncm3::T      | density in n_pcm^3        | | T_K::T           | temperature in K               | | T_eV::T          | temperature in eV              | | P_th_cgs::T      | thermal pressure in Ba         | | P_CR_cgs::T      | cosmic ray pressure in Ba      |\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.GadgetPhysicalUnits","page":"API reference","title":"GadgetUnits.GadgetPhysicalUnits","text":"GadgetPhysicalUnits( DT::DataType, \n                        l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;\n                        a_scale::Real=1.0, hpar::Real=1.0,\n                        γ_th::Real=5.0/3.0, γ_CR::Real=4.0/3.0, xH::Real=0.76)\n\nSet up a unit struct with a given DataType.\n\n\n\n\n\n","category":"type"},{"location":"api/#GadgetUnits.GadgetPhysicalUnits-Union{Tuple{}, Tuple{T}, Tuple{T, T}, Tuple{T, T, T}} where T","page":"API reference","title":"GadgetUnits.GadgetPhysicalUnits","text":"GadgetPhysicalUnits(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;\n                    a_scale::T=1.0, hpar::T=1.0,\n                    γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T\n\nCreates a datatype GadgetPhysicalUnits which holds the conversion factors between comoving code units and physical units. Stores the unit information which can be converted with Unitful or UnitfulAstro.\n\nExamples\n\njulia> GU = GadgetPhysicalUnits()\nGadgetPhysicalUnits(3.085678e21 cm, 100000.0 cm s^-1, 1.989e43 g, 3.085678e16 s, 977.7923542981722 Myr, 1.989e53 erg, 1.2414361549102458e65 eV, 1.0 Gs, 6.769911178294544e-22 g cm^-3, 743179.9340255889 N_e/cm^3, 47.50882854026919 K, 6.769911178294544e-12 Ba, 6.769911178294544e-12 Ba)\njulia> GU.x_cgs\n3.085678e21 cm\n\nKeyword Arguments\n\na_scale::T = 1.0:  Cosmological scale factor of the simulation. Can be passed with the header h as h.time.\nhpar::T = 1.0:     Hubble constant as 'little h'. Can be passed with header h as h.h0.\nγ_th::T = 5.0/3.0: Adiabatic index of gas.\nγ_CR::T = 4.0/3.0: Adiabatic index of cosmic ray component.\nxH::T = 0.76:      Hydrogen fraction of the simulation, if run without chemical model.\n\nFields\n\n| Name                         | Meaning                        | |: –––––––––––––- |:––––––––––––––-  | | x_cgs::Quantity{T}         | position in cm                 | | x_kpc::Quantity{T}         | position in kpc                | | v_cgs::Quantity{T}         | velocity in cm/s               | | v_kms::Quantity{T}         | velocity in km/s               | | m_cgs::Quantity{T}         | mass in g                      | | m_msun::Quantity{T}        | mass in Msun                   | | t_s::Quantity{T}           | time in sec                    | | t_Myr::Quantity{T}         | time in Myr                    | | E_cgs::Quantity{T}         | energy in erg                  | | E_eV::Quantity{T}          | energy in eV                   | | B_cgs::Quantity{T}         | magnetic field in Gauss        | | rho_cgs::Quantity{T}       | density in gcm^3          | | rho_ncm3::Quantity{T}      | density in n_pcm^3        | | T_K::Quantity{T}           | temperature in K               | | T_eV::Quantity{T}          | temperature in eV              | | P_th_cgs::Quantity{T}      | thermal pressure in Ba         | | P_CR_cgs::Quantity{T}      | cosmic ray pressure in Ba      |\n\n\n\n\n\n","category":"method"},{"location":"api/#Private-Functions","page":"API reference","title":"Private Functions","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [GadgetUnits]\nPublic = false\nOrder = [:function]","category":"page"},{"location":"api/#Private-Types","page":"API reference","title":"Private Types","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [GadgetUnits]\nPublic = false\nOrder = [:type]","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"CurrentModule = GadgetUnits\nDocTestSetup = quote\n    using GadgetUnits\nend","category":"page"},{"location":"cosmo/#Cosmological-conversions","page":"Cosmological conversions","title":"Cosmological conversions","text":"","category":"section"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"GadgetUnits.jl supplied some conversion functions to convert observables in given cosmologies. These functions rely on Cosmology.jl, specifically a struct of AbstractCosmology.","category":"page"},{"location":"cosmo/#Cosmology","page":"Cosmological conversions","title":"Cosmology","text":"","category":"section"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"To define an AbstractCosmology from the properties of the simulation you can use the multiple dispatch function","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"c = cosmology(h::SnapshotHeader)","category":"page"},{"location":"cosmo/#Time","page":"Cosmological conversions","title":"Time","text":"","category":"section"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"To get the current age of the universe you can use age","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"t = age(h::SnapshotHeader, units::Bool=true)","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"If the optional parameter units is set to true it will preserve the unit property. Otherwise it will return a Real in units Gyrs.","category":"page"},{"location":"cosmo/#Scale","page":"Cosmological conversions","title":"Scale","text":"","category":"section"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"You can convert arcminutes to kpc for an object at redshift z with armin_to_kpc. This can be used either with a pre-defined Cosmology.AbstractCosmology ","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"arcmin_to_kpc(c::Cosmology.AbstractCosmology, θ::Real, z::Real )","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"or can use a SnapshotHeader struct from GadgetIO.jl which defines the AbstractCosmology on the fly:","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"arcmin_to_kpc(θ::Real, h::SnapshotHeader)","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"If you input θ::Real the function will return kpc, but remove the unit property. Optionally you can define θ::Unitful.AbstractQuantity which preserves the unit property.","category":"page"},{"location":"cosmo/#Radiation","page":"Cosmological conversions","title":"Radiation","text":"","category":"section"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"To convert synchrotron radiation from mJy to W/Hz you can use mJy_to_W. Like in the case of Scale you can use either a pre-defined AbstractCosmology","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"mJy_to_W( c::Cosmology.AbstractCosmology, S::Real, z::Real )","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"or use a SnapshotHeader","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"mJy_to_W(S::Union{Real, Unitful.AbstractQuantity}, h::SnapshotHeader)","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"CurrentModule = GadgetUnits\nDocTestSetup = quote\n    using GadgetUnits\nend","category":"page"},{"location":"conversion_structs/#Conversion-Structs","page":"Basic Unit Conversion","title":"Conversion Structs","text":"","category":"section"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"GadgetUnits.jl uses Unitful.jl and UnitfulAstro.jl to store the unit conversion factors with actual units in place. You can convert the internal units of Gadget into cgs units by defining the object GadgetPhysicalUnits:","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"GU = GadgetPhysicalUnits(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;\n                         a_scale::T=1.0, hpar::T=1.0,\n                         γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"where the keyword arguments are:","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"a_scale::T = 1.0:  Cosmological scale factor of the simulation. Can be passed with the header h as h.time.\nhpar::T = 1.0:     Hubble constant as 'little h'. Can be passed with header h as h.h0.\nγ_th::T = 5.0/3.0: Adiabatic index of gas.\nγ_CR::T = 4.0/3.0: Adiabatic index of cosmic ray component.\nxH::T = 0.76:      Hydrogen fraction of the simulation, if run without chemical model.","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"This returns an object of type GadgetPhysicalUnits.","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"If you want to have the same functionality, but without using Unitful.jl you can construct a similar object GadgetPhysical:","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"GU = GadgetPhysical(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;\n                    a_scale::T=1.0, hpar::T=1.0,\n                    γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"This uses the same conversions, but leaves out the actual unit properties.","category":"page"},{"location":"conversion_structs/#Usage","page":"Basic Unit Conversion","title":"Usage","text":"","category":"section"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"To convert, say positions of gas particles from a cosmological simulation to physical units you can use:","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"\nh     = read_header(filename)\n\npos   = read_snap(filename, \"POS\", 0)\n\nGU    = GadgetPhysicalUnits(a_scale=h.time, hpar=h.h0)\n\npos .*= GU.x_cgs\n","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"If you have different units than the standard Gadget ones you can call the object cunstructor with different values","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"GU = GadgetPhysicalUnits(your_l_unit, your_m_unit, your_v_unit; kwargs...)","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"Converting the units can then be done with Unitful.jl and UnitfulAstro.jl. So if you want to convert the position units from the default cm to Mpc you can do this as:","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"using Unitful\nusing UnitfulAstro\n\npos = read_snap(filename, \"POS\", 0)\npos = @. pos * GU.x_cgs |> u\"Mpc\"","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"If you want to get rid of the units, for example if you need basic datatypes again for a function you can use the funtion ustrip:","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"pos = ustrip(pos)","category":"page"},{"location":"conversion_structs/#Data-Types","page":"Basic Unit Conversion","title":"Data Types","text":"","category":"section"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"You can define the datatype of the struct by passing it as an optional first parameter","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"GadgetPhysicalUnits( DT::DataType, \n                     l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;\n                     a_scale::Real=1.0, hpar::Real=1.0,\n                     γ_th::Real=5.0/3.0, γ_CR::Real=4.0/3.0, xH::Real=0.76)","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"DT defaults to Float64.","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"Please be aware that some unit conversions overflow at Float32!","category":"page"},{"location":"install/#Installing","page":"Install","title":"Installing","text":"","category":"section"},{"location":"install/","page":"Install","title":"Install","text":"As usual with Julia you can install the package via the internal package manager, so via the REPL:","category":"page"},{"location":"install/","page":"Install","title":"Install","text":"julia> ]\npkg> add GadgetUnits","category":"page"},{"location":"install/","page":"Install","title":"Install","text":"Now you should be good to go!","category":"page"},{"location":"","page":"Table of Contents","title":"Table of Contents","text":"CurrentModule = GadgetUnits\nDocTestSetup = quote\n    using GadgetUnits\nend","category":"page"},{"location":"#Table-of-Contents","page":"Table of Contents","title":"Table of Contents","text":"","category":"section"},{"location":"","page":"Table of Contents","title":"Table of Contents","text":"Pages = [ \"index.md\",\n          \"install.md\",\n          \"conversion_structs.md\",\n          \"cosmo.md\",\n          \"api.md\" \n        ]\nDepth = 3","category":"page"}]
}
