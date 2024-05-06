var documenterSearchIndex = {"docs":
[{"location":"api/#API-Reference","page":"API reference","title":"API Reference","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"CurrentModule = GadgetUnits\nDocTestSetup = quote\n    using GadgetUnits\nend","category":"page"},{"location":"api/","page":"API reference","title":"API reference","text":"","category":"page"},{"location":"api/#Exported-Functions","page":"API reference","title":"Exported Functions","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [GadgetUnits]\nPrivate = false\nOrder = [:function]","category":"page"},{"location":"api/#Cosmology.age","page":"API reference","title":"Cosmology.age","text":"age(h::AbstractGadgetHeader, units::Bool=false)\n\nComputes the age of the universe given the properties from AbstractGadgetHeader.\n\n\n\n\n\n","category":"function"},{"location":"api/#Cosmology.cosmology-Tuple{GadgetIO.AbstractGadgetHeader}","page":"API reference","title":"Cosmology.cosmology","text":"cosmology(h::AbstractGadgetHeader)\n\nDefines a Cosmology.AbstractCosmology from the properties of h.\n\n\n\n\n\n","category":"method"},{"location":"api/#Cosmology.lookback_time","page":"API reference","title":"Cosmology.lookback_time","text":"lookback_time(h::AbstractGadgetHeader, units::Bool=false)\n\nComputes the lookback tima given the properties from AbstractGadgetHeader.\n\n\n\n\n\n","category":"function"},{"location":"api/#GadgetUnits.arcmin_to_kpc-Tuple{Cosmology.AbstractCosmology, Real, Real}","page":"API reference","title":"GadgetUnits.arcmin_to_kpc","text":"arcmin_to_kpc(c::Cosmology.AbstractCosmology, θ::Real, z::Real)\n\nConvert arcmin to kpc for a given redshift z and cosmology c. Deletes unit information.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.arcmin_to_kpc-Tuple{Cosmology.AbstractCosmology, Unitful.AbstractQuantity, Real}","page":"API reference","title":"GadgetUnits.arcmin_to_kpc","text":"arcmin_to_kpc(c::Cosmology.AbstractCosmology, θ::Unitful.AbstractQuantity, z::Real)\n\nConvert arcmin to kpc for a given redshift z and cosmology c. Conserves unit information.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.arcmin_to_kpc-Tuple{Real, GadgetIO.AbstractGadgetHeader}","page":"API reference","title":"GadgetUnits.arcmin_to_kpc","text":"arcmin_to_kpc(θ::Real, h::AbstractGadgetHeader)\n\nConvert arcmin to kpc for a given redshift and cosmology taken from AbstractGadgetHeader.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.arcmin_to_kpc-Tuple{Unitful.AbstractQuantity, GadgetIO.AbstractGadgetHeader}","page":"API reference","title":"GadgetUnits.arcmin_to_kpc","text":"arcmin_to_kpc(θ::Unitful.AbstractQuantity, h::AbstractGadgetHeader)\n\nConvert arcmin to kpc for a given redshift and cosmology taken from AbstractGadgetHeader.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.mJy_to_W-Tuple{Cosmology.AbstractCosmology, Real, Real}","page":"API reference","title":"GadgetUnits.mJy_to_W","text":"mJy_to_W( c::Cosmology.AbstractCosmology, S::Real, z::Real )\n\nConverts synchrotron flux density S in [mJy] to [W/Hz] for a given redshift z and cosmology c. Deletes unit information.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.mJy_to_W-Tuple{Cosmology.AbstractCosmology, Unitful.AbstractQuantity, Real}","page":"API reference","title":"GadgetUnits.mJy_to_W","text":"mJy_to_W( c::Cosmology.AbstractCosmology, S::Unitful.AbstractQuantity, z::Real )\n\nConverts synchrotron flux density S in [mJy] to [W/Hz] for a given redshift z and cosmology c. Conserves unit information.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.mJy_to_W-Tuple{Real, Real}","page":"API reference","title":"GadgetUnits.mJy_to_W","text":"mJy_to_W(S::Union{Real,Unitful.AbstractQuantity}, d::Real)\n\nConverts synchrotron flux density S in [mJy] to [W/Hz] for a given distance d in [Mpc].\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.mJy_to_W-Tuple{Union{Real, Unitful.AbstractQuantity}, GadgetIO.AbstractGadgetHeader}","page":"API reference","title":"GadgetUnits.mJy_to_W","text":"mJy_to_W(S::Union{Real, Unitful.AbstractQuantity}, h::AbstractGadgetHeader)\n\nConverts synchrotron flux density S in [mJy] to [W/Hz] for a given redshift and cosmology taken from AbstractGadgetHeader.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.mJy_to_W-Tuple{Unitful.AbstractQuantity, Real}","page":"API reference","title":"GadgetUnits.mJy_to_W","text":"mJy_to_W(S::Union{Real,Unitful.AbstractQuantity}, d::Real)\n\nConverts synchrotron flux density S in [mJy] to [W/Hz] for a given distance d in [Mpc].\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.mJy_to_W-Tuple{Unitful.AbstractQuantity, Unitful.AbstractQuantity}","page":"API reference","title":"GadgetUnits.mJy_to_W","text":"mJy_to_W(S::Union{Real,Unitful.AbstractQuantity}, d::Unitful.AbstractQuantity)\n\nConverts synchrotron flux density S in [mJy] to [W/Hz] for a given distance d.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.redshift-Tuple{Real, Cosmology.AbstractCosmology}","page":"API reference","title":"GadgetUnits.redshift","text":"redshift(c::AbstractCosmology, t_Gyrs::Real)\n\nComputes the redshift at a given age of the universe in Gyrs for a cosmology t.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.redshift-Tuple{Real, GadgetIO.AbstractGadgetHeader}","page":"API reference","title":"GadgetUnits.redshift","text":"redshift(t::Real, h::AbstractGadgetHeader)\n\nSee redshift.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.redshift-Tuple{Unitful.AbstractQuantity, Cosmology.AbstractCosmology}","page":"API reference","title":"GadgetUnits.redshift","text":"redshift(t::Unitful.AbstractQuantity, h::AbstractGadgetHeader)\n\nSee redshift.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.redshift-Tuple{Unitful.AbstractQuantity, GadgetIO.AbstractGadgetHeader}","page":"API reference","title":"GadgetUnits.redshift","text":"redshift(t::Unitful.AbstractQuantity, h::AbstractGadgetHeader)\n\nSee redshift.\n\n\n\n\n\n","category":"method"},{"location":"api/#Exported-Types","page":"API reference","title":"Exported Types","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [GadgetUnits]\nPrivate = false\nOrder = [:type]","category":"page"},{"location":"api/#GadgetUnits.GadgetPhysical","page":"API reference","title":"GadgetUnits.GadgetPhysical","text":"GadgetPhysical( DT::DataType, \n                l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;\n                a_scale::Real=1.0, hpar::Real=1.0,\n                γ_th::Real=5.0/3.0, xH::Real=0.752)\n\nSet up a unit struct with a given DataType.\n\n\n\n\n\n","category":"type"},{"location":"api/#GadgetUnits.GadgetPhysical-Union{Tuple{GadgetIO.SnapshotHeader}, Tuple{T}, Tuple{GadgetIO.SnapshotHeader, T}, Tuple{GadgetIO.SnapshotHeader, T, T}, Tuple{GadgetIO.SnapshotHeader, T, T, T}} where T","page":"API reference","title":"GadgetUnits.GadgetPhysical","text":"GadgetPhysical( h::SnapshotHeader, \n                l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;\n                γ_th::Real=5.0/3.0, xH::Real=0.752)\n\nSet up a unit struct with a given SnapshotHeader. Only works in 64-bits.\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.GadgetPhysical-Union{Tuple{}, Tuple{T}, Tuple{T, T}, Tuple{T, T, T}} where T","page":"API reference","title":"GadgetUnits.GadgetPhysical","text":"GadgetPhysical(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;\n               a_scale::T=1.0, hpar::T=1.0,\n               γ_th::T=5.0/3.0, xH::T=0.752) where T\n\nCreates a struct GadgetPhysical which holds the conversion factors between comoving code units and physical units, without unit information.\n\nKeyword arugments specify:\n\nArguments\n\na_scale::T = 1.0:  Cosmological scale factor of the simulation. Can be passed with the header h as h.time.\nhpar::T = 1.0:     Hubble constant as 'little h'. Can be passed with header h as h.h0.\nγ_th::T = 5.0/3.0: Adiabatic index of gas.\nxH::T = 0.752:      Hydrogen fraction of the simulation, if run without chemical model.\n\nFields\n\nName Meaning\nx_cgs::T position in cm\nx_physical::T position in kpc\nv_cgs::T velocity in cm/s\nv_physical::T velocity in km/s\nm_cgs::T mass in g\nm_msun::T mass in Msun\nm_physical::T mass in 10^10 Msun\nt_s::T time in sec\nt_Myr::T time in Myr\nE_cgs::T energy in erg\nE_eV::T energy in eV\nB_cgs::T magnetic field in Gauss\nrho_physical::T density in 10^10 Msun/kpc^3\nrho_cgs::T density in gcm^3\nrho_ncm3::T density in n_pcm^3\nT_K::T temperature in K\nT_eV::T temperature in eV\nP_th_cgs::T thermal pressure in Ba\nP_CR_cgs::T cosmic ray pressure in Ba\nCR_norm::T LMBSPECTRALCRs Norm in cgs units\n\n\n\n\n\n","category":"method"},{"location":"api/#GadgetUnits.GadgetPhysicalUnits","page":"API reference","title":"GadgetUnits.GadgetPhysicalUnits","text":"GadgetPhysicalUnits( DT::DataType, \n                        l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;\n                        a_scale::Real=1.0, hpar::Real=1.0,\n                        γ_th::Real=5.0/3.0, xH::Real=0.752)\n\nSet up a unit struct with a given DataType.\n\n\n\n\n\n","category":"type"},{"location":"api/#GadgetUnits.GadgetPhysicalUnits-2","page":"API reference","title":"GadgetUnits.GadgetPhysicalUnits","text":"GadgetPhysical( h::SnapshotHeader, \n                l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;\n                γ_th::Real=5.0/3.0, xH::Real=0.752)\n\nSet up a unit struct with a given SnapshotHeader. Only works in 64-bits.\n\n\n\n\n\n","category":"type"},{"location":"api/#GadgetUnits.GadgetPhysicalUnits-Union{Tuple{}, Tuple{T}, Tuple{T, T}, Tuple{T, T, T}} where T","page":"API reference","title":"GadgetUnits.GadgetPhysicalUnits","text":"GadgetPhysicalUnits(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;\n                    a_scale::T=1.0, hpar::T=1.0,\n                    γ_th::T=5.0/3.0, xH::T=0.752) where T\n\nCreates a struct GadgetPhysicalUnits which holds the conversion factors between comoving code units and physical units. Stores the unit information which can be converted with Unitful or UnitfulAstro.\n\nKeyword Arguments\n\na_scale::T = 1.0:  Cosmological scale factor of the simulation. Can be passed with the header h as h.time.\nhpar::T = 1.0:     Hubble constant as 'little h'. Can be passed with header h as h.h0.\nγ_th::T = 5.0/3.0: Adiabatic index of gas.\nxH::T = 0.752:      Hydrogen fraction of the simulation, if run without chemical model.\n\nFields\n\nName Meaning\nx_cgs::Quantity{T} position in cm\nx_physical::Quantity{T} position in kpc\nv_cgs::Quantity{T} velocity in cm/s\nv_physical::Quantity{T} velocity in km/s\nm_cgs::Quantity{T} mass in g\nm_msun::Quantity{T} mass in solar masses\nt_s::Quantity{T} time in sec\nt_Myr::Quantity{T} time in Myr\nE_cgs::Quantity{T} energy in erg\nE_eV::Quantity{T} energy in eV\nB_cgs::Quantity{T} magnetic field in Gauss\nrho_cgs::Quantity{T} density in gcm^3\nrho_ncm3::Quantity{T} density in n_pcm^3\nT_K::Quantity{T} temperature in K\nT_eV::Quantity{T} temperature in eV\nP_th_cgs::Quantity{T} thermal pressure in Ba\nP_CR_cgs::Quantity{T} cosmic ray pressure in Ba\nCR_norm::Quantity{T} LMBSPECTRALCRs Norm in cgs units\n\n\n\n\n\n","category":"method"},{"location":"api/#Private-Functions","page":"API reference","title":"Private Functions","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [GadgetUnits]\nPublic = false\nOrder = [:function]","category":"page"},{"location":"api/#Private-Types","page":"API reference","title":"Private Types","text":"","category":"section"},{"location":"api/","page":"API reference","title":"API reference","text":"Modules = [GadgetUnits]\nPublic = false\nOrder = [:type]","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"CurrentModule = GadgetUnits\nDocTestSetup = quote\n    using GadgetUnits\nend","category":"page"},{"location":"cosmo/#Cosmological-conversions","page":"Cosmological conversions","title":"Cosmological conversions","text":"","category":"section"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"GadgetUnits.jl supplied some conversion functions to convert observables in given cosmologies. These functions rely on Cosmology.jl, specifically a struct of AbstractCosmology.","category":"page"},{"location":"cosmo/#Cosmology","page":"Cosmological conversions","title":"Cosmology","text":"","category":"section"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"To define an AbstractCosmology from the properties of the simulation you can use the multiple dispatch function","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"cosmology","category":"page"},{"location":"cosmo/#Cosmology.cosmology","page":"Cosmological conversions","title":"Cosmology.cosmology","text":"cosmology(h::AbstractGadgetHeader)\n\nDefines a Cosmology.AbstractCosmology from the properties of h.\n\n\n\n\n\n","category":"function"},{"location":"cosmo/#Time","page":"Cosmological conversions","title":"Time","text":"","category":"section"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"To get the current age of the universe you can use age","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"age","category":"page"},{"location":"cosmo/#Cosmology.age","page":"Cosmological conversions","title":"Cosmology.age","text":"age(h::AbstractGadgetHeader, units::Bool=false)\n\nComputes the age of the universe given the properties from AbstractGadgetHeader.\n\n\n\n\n\n","category":"function"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"If the optional parameter units is set to true it will preserve the unit property. Otherwise it will return a Real in units Gyrs.","category":"page"},{"location":"cosmo/#Scale","page":"Cosmological conversions","title":"Scale","text":"","category":"section"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"You can convert arcminutes to kpc for an object at redshift z with arcmin_to_kpc. This can be used either with a pre-defined Cosmology.AbstractCosmology or can use a AbstractGadgetHeader struct from GadgetIO.jl which defines the AbstractCosmology on the fly. If you input θ::Real the function will return kpc, but remove the unit property. Optionally you can define θ::Unitful.AbstractQuantity which preserves the unit property.","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"arcmin_to_kpc","category":"page"},{"location":"cosmo/#GadgetUnits.arcmin_to_kpc","page":"Cosmological conversions","title":"GadgetUnits.arcmin_to_kpc","text":"arcmin_to_kpc(c::Cosmology.AbstractCosmology, θ::Unitful.AbstractQuantity, z::Real)\n\nConvert arcmin to kpc for a given redshift z and cosmology c. Conserves unit information.\n\n\n\n\n\narcmin_to_kpc(c::Cosmology.AbstractCosmology, θ::Real, z::Real)\n\nConvert arcmin to kpc for a given redshift z and cosmology c. Deletes unit information.\n\n\n\n\n\narcmin_to_kpc(θ::Real, h::AbstractGadgetHeader)\n\nConvert arcmin to kpc for a given redshift and cosmology taken from AbstractGadgetHeader.\n\n\n\n\n\narcmin_to_kpc(θ::Unitful.AbstractQuantity, h::AbstractGadgetHeader)\n\nConvert arcmin to kpc for a given redshift and cosmology taken from AbstractGadgetHeader.\n\n\n\n\n\n","category":"function"},{"location":"cosmo/#Radiation","page":"Cosmological conversions","title":"Radiation","text":"","category":"section"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"To convert synchrotron radiation from mJy to W/Hz you can use mJy_to_W. Like in the case of Scale you can use either a pre-defined AbstractCosmology or use an AbstractGadgetHeader.","category":"page"},{"location":"cosmo/","page":"Cosmological conversions","title":"Cosmological conversions","text":"mJy_to_W","category":"page"},{"location":"cosmo/#GadgetUnits.mJy_to_W","page":"Cosmological conversions","title":"GadgetUnits.mJy_to_W","text":"mJy_to_W( c::Cosmology.AbstractCosmology, S::Unitful.AbstractQuantity, z::Real )\n\nConverts synchrotron flux density S in [mJy] to [W/Hz] for a given redshift z and cosmology c. Conserves unit information.\n\n\n\n\n\nmJy_to_W( c::Cosmology.AbstractCosmology, S::Real, z::Real )\n\nConverts synchrotron flux density S in [mJy] to [W/Hz] for a given redshift z and cosmology c. Deletes unit information.\n\n\n\n\n\nmJy_to_W(S::Union{Real, Unitful.AbstractQuantity}, h::AbstractGadgetHeader)\n\nConverts synchrotron flux density S in [mJy] to [W/Hz] for a given redshift and cosmology taken from AbstractGadgetHeader.\n\n\n\n\n\nmJy_to_W(S::Union{Real,Unitful.AbstractQuantity}, d::Unitful.AbstractQuantity)\n\nConverts synchrotron flux density S in [mJy] to [W/Hz] for a given distance d.\n\n\n\n\n\nmJy_to_W(S::Union{Real,Unitful.AbstractQuantity}, d::Real)\n\nConverts synchrotron flux density S in [mJy] to [W/Hz] for a given distance d in [Mpc].\n\n\n\n\n\nmJy_to_W(S::Union{Real,Unitful.AbstractQuantity}, d::Real)\n\nConverts synchrotron flux density S in [mJy] to [W/Hz] for a given distance d in [Mpc].\n\n\n\n\n\n","category":"function"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"CurrentModule = GadgetUnits\nDocTestSetup = quote\n    using GadgetUnits\nend","category":"page"},{"location":"conversion_structs/#Conversion-Structs","page":"Basic Unit Conversion","title":"Conversion Structs","text":"","category":"section"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"GadgetUnits.jl uses Unitful.jl and UnitfulAstro.jl to store the unit conversion factors with actual units in place. You can convert the internal units of Gadget into cgs units by defining the struct","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"GadgetPhysicalUnits","category":"page"},{"location":"conversion_structs/#GadgetUnits.GadgetPhysicalUnits","page":"Basic Unit Conversion","title":"GadgetUnits.GadgetPhysicalUnits","text":"GadgetPhysicalUnits(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;\n                    a_scale::T=1.0, hpar::T=1.0,\n                    γ_th::T=5.0/3.0, xH::T=0.752) where T\n\nCreates a struct GadgetPhysicalUnits which holds the conversion factors between comoving code units and physical units. Stores the unit information which can be converted with Unitful or UnitfulAstro.\n\nKeyword Arguments\n\na_scale::T = 1.0:  Cosmological scale factor of the simulation. Can be passed with the header h as h.time.\nhpar::T = 1.0:     Hubble constant as 'little h'. Can be passed with header h as h.h0.\nγ_th::T = 5.0/3.0: Adiabatic index of gas.\nxH::T = 0.752:      Hydrogen fraction of the simulation, if run without chemical model.\n\nFields\n\nName Meaning\nx_cgs::Quantity{T} position in cm\nx_physical::Quantity{T} position in kpc\nv_cgs::Quantity{T} velocity in cm/s\nv_physical::Quantity{T} velocity in km/s\nm_cgs::Quantity{T} mass in g\nm_msun::Quantity{T} mass in solar masses\nt_s::Quantity{T} time in sec\nt_Myr::Quantity{T} time in Myr\nE_cgs::Quantity{T} energy in erg\nE_eV::Quantity{T} energy in eV\nB_cgs::Quantity{T} magnetic field in Gauss\nrho_cgs::Quantity{T} density in gcm^3\nrho_ncm3::Quantity{T} density in n_pcm^3\nT_K::Quantity{T} temperature in K\nT_eV::Quantity{T} temperature in eV\nP_th_cgs::Quantity{T} thermal pressure in Ba\nP_CR_cgs::Quantity{T} cosmic ray pressure in Ba\nCR_norm::Quantity{T} LMBSPECTRALCRs Norm in cgs units\n\n\n\n\n\nGadgetPhysicalUnits( DT::DataType, \n                        l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;\n                        a_scale::Real=1.0, hpar::Real=1.0,\n                        γ_th::Real=5.0/3.0, xH::Real=0.752)\n\nSet up a unit struct with a given DataType.\n\n\n\n\n\nGadgetPhysical( h::SnapshotHeader, \n                l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;\n                γ_th::Real=5.0/3.0, xH::Real=0.752)\n\nSet up a unit struct with a given SnapshotHeader. Only works in 64-bits.\n\n\n\n\n\n","category":"type"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"If you want to have the same functionality, but without using Unitful.jl you can construct a similar struct:","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"GadgetPhysical","category":"page"},{"location":"conversion_structs/#GadgetUnits.GadgetPhysical","page":"Basic Unit Conversion","title":"GadgetUnits.GadgetPhysical","text":"GadgetPhysical(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;\n               a_scale::T=1.0, hpar::T=1.0,\n               γ_th::T=5.0/3.0, xH::T=0.752) where T\n\nCreates a struct GadgetPhysical which holds the conversion factors between comoving code units and physical units, without unit information.\n\nKeyword arugments specify:\n\nArguments\n\na_scale::T = 1.0:  Cosmological scale factor of the simulation. Can be passed with the header h as h.time.\nhpar::T = 1.0:     Hubble constant as 'little h'. Can be passed with header h as h.h0.\nγ_th::T = 5.0/3.0: Adiabatic index of gas.\nxH::T = 0.752:      Hydrogen fraction of the simulation, if run without chemical model.\n\nFields\n\nName Meaning\nx_cgs::T position in cm\nx_physical::T position in kpc\nv_cgs::T velocity in cm/s\nv_physical::T velocity in km/s\nm_cgs::T mass in g\nm_msun::T mass in Msun\nm_physical::T mass in 10^10 Msun\nt_s::T time in sec\nt_Myr::T time in Myr\nE_cgs::T energy in erg\nE_eV::T energy in eV\nB_cgs::T magnetic field in Gauss\nrho_physical::T density in 10^10 Msun/kpc^3\nrho_cgs::T density in gcm^3\nrho_ncm3::T density in n_pcm^3\nT_K::T temperature in K\nT_eV::T temperature in eV\nP_th_cgs::T thermal pressure in Ba\nP_CR_cgs::T cosmic ray pressure in Ba\nCR_norm::T LMBSPECTRALCRs Norm in cgs units\n\n\n\n\n\nGadgetPhysical( DT::DataType, \n                l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;\n                a_scale::Real=1.0, hpar::Real=1.0,\n                γ_th::Real=5.0/3.0, xH::Real=0.752)\n\nSet up a unit struct with a given DataType.\n\n\n\n\n\nGadgetPhysical( h::SnapshotHeader, \n                l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;\n                γ_th::Real=5.0/3.0, xH::Real=0.752)\n\nSet up a unit struct with a given SnapshotHeader. Only works in 64-bits.\n\n\n\n\n\n","category":"type"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"This uses the same conversions, but leaves out the actual unit properties.","category":"page"},{"location":"conversion_structs/#Usage","page":"Basic Unit Conversion","title":"Usage","text":"","category":"section"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"To convert, say positions of gas particles from a cosmological simulation to physical units you can use:","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"h     = read_header(filename)\npos   = read_snap(filename, \"POS\", 0)\nGU    = GadgetPhysicalUnits(a_scale=h.time, hpar=h.h0)\npos .*= GU.x_cgs","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"If you have different units than the standard Gadget ones you can call the cunstructor with different values","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"GU = GadgetPhysicalUnits(your_l_unit, your_m_unit, your_v_unit; kwargs...)","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"Converting the units can then be done with Unitful.jl and UnitfulAstro.jl. So if you want to convert the position units from the default cm to Mpc you can do this as:","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"using Unitful\nusing UnitfulAstro\n\npos = read_snap(filename, \"POS\", 0)\npos = @. pos * GU.x_cgs |> u\"Mpc\"","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"If you want to get rid of the units, for example if you need basic datatypes again for a function you can use the funtion ustrip:","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"pos = ustrip(pos)","category":"page"},{"location":"conversion_structs/#Data-Types","page":"Basic Unit Conversion","title":"Data Types","text":"","category":"section"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"You can define the datatype of the struct by passing it as an optional first parameter","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"GadgetPhysicalUnits(::DataType)","category":"page"},{"location":"conversion_structs/#GadgetUnits.GadgetPhysicalUnits-Tuple{DataType}","page":"Basic Unit Conversion","title":"GadgetUnits.GadgetPhysicalUnits","text":"GadgetPhysicalUnits(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;\n                    a_scale::T=1.0, hpar::T=1.0,\n                    γ_th::T=5.0/3.0, xH::T=0.752) where T\n\nCreates a struct GadgetPhysicalUnits which holds the conversion factors between comoving code units and physical units. Stores the unit information which can be converted with Unitful or UnitfulAstro.\n\nKeyword Arguments\n\na_scale::T = 1.0:  Cosmological scale factor of the simulation. Can be passed with the header h as h.time.\nhpar::T = 1.0:     Hubble constant as 'little h'. Can be passed with header h as h.h0.\nγ_th::T = 5.0/3.0: Adiabatic index of gas.\nxH::T = 0.752:      Hydrogen fraction of the simulation, if run without chemical model.\n\nFields\n\nName Meaning\nx_cgs::Quantity{T} position in cm\nx_physical::Quantity{T} position in kpc\nv_cgs::Quantity{T} velocity in cm/s\nv_physical::Quantity{T} velocity in km/s\nm_cgs::Quantity{T} mass in g\nm_msun::Quantity{T} mass in solar masses\nt_s::Quantity{T} time in sec\nt_Myr::Quantity{T} time in Myr\nE_cgs::Quantity{T} energy in erg\nE_eV::Quantity{T} energy in eV\nB_cgs::Quantity{T} magnetic field in Gauss\nrho_cgs::Quantity{T} density in gcm^3\nrho_ncm3::Quantity{T} density in n_pcm^3\nT_K::Quantity{T} temperature in K\nT_eV::Quantity{T} temperature in eV\nP_th_cgs::Quantity{T} thermal pressure in Ba\nP_CR_cgs::Quantity{T} cosmic ray pressure in Ba\nCR_norm::Quantity{T} LMBSPECTRALCRs Norm in cgs units\n\n\n\n\n\nGadgetPhysicalUnits( DT::DataType, \n                        l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;\n                        a_scale::Real=1.0, hpar::Real=1.0,\n                        γ_th::Real=5.0/3.0, xH::Real=0.752)\n\nSet up a unit struct with a given DataType.\n\n\n\n\n\n","category":"method"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"DT defaults to Float64.","category":"page"},{"location":"conversion_structs/","page":"Basic Unit Conversion","title":"Basic Unit Conversion","text":"Please be aware that some unit conversions overflow at Float32!","category":"page"},{"location":"install/#Installing","page":"Install","title":"Installing","text":"","category":"section"},{"location":"install/","page":"Install","title":"Install","text":"As usual with Julia you can install the package via the internal package manager, so via the REPL:","category":"page"},{"location":"install/","page":"Install","title":"Install","text":"julia> ]\npkg> add GadgetUnits","category":"page"},{"location":"install/","page":"Install","title":"Install","text":"Now you should be good to go!","category":"page"},{"location":"","page":"Table of Contents","title":"Table of Contents","text":"CurrentModule = GadgetUnits\nDocTestSetup = quote\n    using GadgetUnits\nend","category":"page"},{"location":"#GadgetUnits.jl","page":"Table of Contents","title":"GadgetUnits.jl","text":"","category":"section"},{"location":"","page":"Table of Contents","title":"Table of Contents","text":"This package provides some basic unit conversion functionality to work with the SPH code Gadget2 by Volker Springel. It is based on the unit conversions listed by Klaus Dolag on his HowTo page (restricted). Further functionality provides tools for cosmological simulations.","category":"page"},{"location":"","page":"Table of Contents","title":"Table of Contents","text":"If you use GadgetUnits.jl in publications please cite it.","category":"page"},{"location":"#Table-of-Contents","page":"Table of Contents","title":"Table of Contents","text":"","category":"section"},{"location":"","page":"Table of Contents","title":"Table of Contents","text":"Pages = [ \"index.md\",\n          \"install.md\",\n          \"conversion_structs.md\",\n          \"cosmo.md\",\n          \"api.md\" \n        ]\nDepth = 3","category":"page"}]
}
