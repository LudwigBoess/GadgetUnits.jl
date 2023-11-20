
struct GadgetPhysicalUnits{T}

    x_cgs::Quantity{T}        # position in cm
    x_physical::Quantity{T}        # position in kpc

    v_cgs::Quantity{T}     # velocity in cm/s
    v_physical::Quantity{T}     # velocity in km/s

    m_cgs::Quantity{T}     # mass in g
    m_msun::Quantity{T}     # mass in Msun

    t_s::Quantity{T}          # time in sec
    t_Myr::Quantity{T}        # time in Myr

    E_cgs::Quantity{T}        # energy in erg
    E_eV::Quantity{T}         # energy in eV

    B_cgs::Quantity{T}     # magnetic field in Gauss

    rho_physical::Quantity{T}   # density in 10^10 Msun/kpc^3
    rho_cgs::Quantity{T}   # density in g/cm^3
    rho_ncm3::Quantity{T}    # density in mp/cm^3

    T_K::Quantity{T}           # temperature in K
    T_eV::Quantity{T}            # temperature in eV

    P_th_cgs::Quantity{T}      # thermal pressure in Ba
    P_CR_cgs::Quantity{T}      # cosmic ray pressure in Ba

    CR_norm::Quantity{T}       # LMB_SPECTRAL_CRs Norm in cgs units
end


"""
    GadgetPhysicalUnits(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;
                        a_scale::T=1.0, hpar::T=1.0,
                        γ_th::T=5.0/3.0, xH::T=0.752) where T

Creates a struct `GadgetPhysicalUnits` which holds the conversion factors between comoving code units and physical units.
Stores the unit information which can be converted with Unitful or UnitfulAstro.

# Keyword Arguments
- `a_scale::T = 1.0`:  Cosmological scale factor of the simulation. Can be passed with the header `h` as `h.time`.
- `hpar::T = 1.0`:     Hubble constant as 'little h'. Can be passed with header `h` as `h.h0`.
- `γ_th::T = 5.0/3.0`: Adiabatic index of gas.
- `xH::T = 0.752`:      Hydrogen fraction of the simulation, if run without chemical model.

# Fields
| Name                         | Meaning                        |
|:---------------------------- |:-----------------------------  |
| `x_cgs::Quantity{T}`         | position in cm                 |
| `x_physical::Quantity{T}`    | position in kpc                |
| `v_cgs::Quantity{T}`         | velocity in cm/s               |
| `v_physical::Quantity{T}`    | velocity in km/s               |
| `m_cgs::Quantity{T}`         | mass in g                      |
| `m_msun::Quantity{T}`        | mass in solar masses           |
| `t_s::Quantity{T}`           | time in sec                    |
| `t_Myr::Quantity{T}`         | time in Myr                    |
| `E_cgs::Quantity{T}`         | energy in erg                  |
| `E_eV::Quantity{T}`          | energy in eV                   |
| `B_cgs::Quantity{T}`         | magnetic field in Gauss        |
| `rho_cgs::Quantity{T}`       | density in ``g/cm^3``          |
| `rho_ncm3::Quantity{T}`      | density in ``n_p/cm^3``        |
| `T_K::Quantity{T}`           | temperature in K               |
| `T_eV::Quantity{T}`          | temperature in eV              |
| `P_th_cgs::Quantity{T}`      | thermal pressure in Ba         |
| `P_CR_cgs::Quantity{T}`      | cosmic ray pressure in Ba      |
| `CR_norm::Quantity{T}`       | LMB_SPECTRAL_CRs Norm in cgs units |
"""
function GadgetPhysicalUnits(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;
                             a_scale::T=1.0, hpar::T=1.0,
                             γ_th::T=5.0 / 3.0, xH::T=0.752) where {T}


    GU = GadgetPhysical(l_unit, m_unit, v_unit, a_scale=a_scale, hpar=hpar,
                        γ_th=γ_th, xH=xH)

    # convert comoving output to physical units
    x_cgs      = GU.x_cgs * 1.0u"cm"
    x_physical = GU.x_physical * 1.0u"kpc"

    v_cgs      = GU.v_cgs * 1.0u"cm/s"
    v_physical = GU.v_physical * 1.0u"km/s"

    m_cgs      = GU.m_cgs  * 1.0u"g"
    m_msun     = GU.m_msun * 1.0u"Msun"

    t_s        = GU.t_s * 1.0u"s"
    t_Myr      = t_s |> u"Myr"

    E_cgs = GU.E_cgs * 1.0u"erg"
    E_eV = E_cgs |> u"eV"

    B_cgs = 1.0u"g^(1/2) / cm^(1/2) / s"    # gadget outputs in Gauss

    rho_cgs = GU.rho_cgs * 1.0u"g/cm^3"
    rho_physical = rho_cgs |> u"Msun/kpc^3"

    rho_ncm3 = GU.rho_ncm3 * 1.0u"cm^-3"

    T_cgs = GU.T_K * 1.0u"K"
    T_eV  = GU.T_eV * 1.0u"eV"

    P_th_cgs = GU.P_th_cgs * 1.0u"erg/cm^3"
    P_CR_cgs = GU.P_CR_cgs * 1.0u"erg/cm^3"

    CR_norm = GU.CR_norm * 1.0u"erg/cm^3"

    GadgetPhysicalUnits{T}(x_cgs, x_physical,
        v_cgs, v_physical,
        m_cgs, m_msun,
        t_s, t_Myr,
        E_cgs, E_eV,
        B_cgs,
        rho_physical, rho_cgs, rho_ncm3,
        T_cgs, T_eV,
        P_th_cgs,
        P_CR_cgs,
        CR_norm
        )

end

"""
    GadgetPhysicalUnits( DT::DataType, 
                            l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;
                            a_scale::Real=1.0, hpar::Real=1.0,
                            γ_th::Real=5.0/3.0, xH::Real=0.752)

Set up a unit struct with a given `DataType`.
"""
GadgetPhysicalUnits( DT::DataType, 
                        l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;
                        a_scale::Real=1.0, hpar::Real=1.0,
                        γ_th::Real=5.0/3.0, xH::Real=0.752) = 
    GadgetPhysicalUnits(DT(l_unit), DT(m_unit), DT(v_unit),
                        a_scale=DT(a_scale), hpar=DT(hpar), 
                        γ_th=DT(γ_th), xH=DT(xH))



"""
    GadgetPhysical( h::SnapshotHeader, 
                    l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;
                    γ_th::Real=5.0/3.0, xH::Real=0.752)

Set up a unit struct with a given `SnapshotHeader`.
Only works in 64-bits.
"""
function GadgetPhysicalUnits( h::SnapshotHeader, 
                            l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;
                            γ_th::Real=5.0/3.0, xH::Real=0.752)

    # calculate scale factor from header -> if non-cosmo sim: h.z = 0.0 for every snapshot    
    a_scale = 1 / ( 1 + h.z )

    GadgetPhysicalUnits( l_unit, m_unit, v_unit,
                a_scale=a_scale, hpar=h.h0, 
                γ_th=γ_th, xH=xH )
end