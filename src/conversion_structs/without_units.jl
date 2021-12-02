struct GadgetPhysical{T}

    x_cgs::T         # position in cm
    x_physical::T    # position in kpc

    v_cgs::T         # velocity in cm/s
    v_physical::T    # velocity in km/s

    m_cgs::T         # mass in g
    m_msun::T        # mass in Msun
    m_physical::T    # mass in 10^10 Msun

    t_s::T           # time in sec
    t_Myr::T         # time in Myr

    E_cgs::T         # energy in erg
    E_eV::T          # energy in eV

    B_cgs::T         # magnetic field in Gauss

    rho_physical::T   # density in 10^10 Msun/kpc^3
    rho_cgs::T       # density in g/cm^3
    rho_ncm3::T      # density in N_p/cm^3

    T_K::T           # temperature in K
    T_eV::T          # temperature in eV

    P_th_cgs::T      # thermal pressure in Ba
    P_CR_cgs::T      # cosmic ray pressure in Ba
end

"""
    GadgetPhysical(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;
                a_scale::T=1.0, hpar::T=1.0,
                γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T

Creates a datatype GadgetPhysical which holds the conversion factors between comoving code units and physical units, without unit information.

Keyword arugments specify:
# Arguments
- `a_scale::T = 1.0`:  Cosmological scale factor of the simulation. Can be passed with the header `h` as `h.time`.
- `hpar::T = 1.0`:     Hubble constant as 'little h'. Can be passed with header `h` as `h.h0`.
- `γ_th::T = 5.0/3.0`: Adiabatic index of gas.
- `γ_CR::T = 4.0/3.0`: Adiabatic index of cosmic ray component.
- `xH::T = 0.76`:      Hydrogen fraction of the simulation, if run without chemical model.

# Fields

|Name | Meaning |
|:--- |:--- |
|`x_cgs::T` | position in cm |
|`x_physical::T` | position in kpc |
|`v_cgs::T` | velocity in cm/s |
|`v_physical::T` | velocity in km/s |
|`m_cgs::T` | mass in g |
|`m_msun::T` | mass in Msun |
|`m_physical::T` | mass in 10^10 Msun |
|`t_s::T` | time in sec |
|`t_Myr::T` | time in Myr |
|`E_cgs::T` | energy in erg |
|`E_eV::T` | energy in eV |
|`B_cgs::T` | magnetic field in Gauss |
|`rho_physical::T` | density in 10^10 Msun/kpc^3 |
|`rho_cgs::T` | density in ``g/cm^3`` |
|`rho_ncm3::T` | density in ``n_p/cm^3`` |
|`T_K::T` | temperature in K |
|`T_eV::T` | temperature in eV |
|`P_th_cgs::T` | thermal pressure in Ba |
|`P_CR_cgs::T` | cosmic ray pressure in Ba |

"""
function GadgetPhysical(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;
                        a_scale::T=1.0, hpar::T=1.0,
                        γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T

    # some basic constants
    kB      = 1.38066e-16
    mp      = 1.6726e-24
    eV2cgs  = 1.602e-12

    yhelium = ( 1.0 - xH ) / ( 4.0 * xH )
    mean_mol_weight = (1.0 + 4.0 * yhelium) / (1.0 + 3.0 * yhelium + 1.0)

    n2ne =  (  xH + 0.5 * ( 1.0 - xH ) ) / 
            ( 2.0 * xH + 0.75 * ( 1.0 - xH ) )       # conversion n_pat -> n_electrons
    umu  =  4.0 / ( 5.0 * xH + 3.0 )                 # mean molucular weight in hydr. mass

    # convert comoving output to physical units
    x_cgs      = l_unit * a_scale / hpar
    x_physical = a_scale / hpar
    
    v_cgs      = v_unit * sqrt(a_scale)
    v_physical = sqrt(a_scale)

    m_cgs      = m_unit / hpar
    m_physical = 1.0 / hpar

    t_unit  = l_unit / v_unit
    t_s     = t_unit * sqrt(a_scale) / hpar  # in sec
    t_Myr   = t_s / 3.15576e13

    E_cgs = m_cgs * v_cgs^2
    E_eV  = E_cgs * 6.242e+11

    B_cgs = 1.0   # gadget outputs in cgs

    rho_physical = m_physical / x_physical^3
    rho_cgs      = m_unit/l_unit^3 * hpar^2 / a_scale^3
    rho_ncm3     = rho_cgs * n2ne/( umu * mp )


    T_cgs = (γ_th - 1.0) * v_unit^2 * mean_mol_weight * mp / kB
    T_eV  = T_cgs * kB / eV2cgs

    P_th_cgs = a_scale^(-3) * E_cgs / l_unit^3 * hpar^2
    P_CR_cgs = a_scale^(-4) * E_cgs / l_unit^3 * hpar^2

    GadgetPhysical{T}(x_cgs, x_physical,
        v_cgs, v_physical,
        m_cgs, m_cgs/T(1.989e33), m_physical,
        t_s, t_Myr,
        E_cgs, E_eV,
        B_cgs,
        rho_physical, rho_cgs, rho_ncm3,
        T_cgs, T_eV,
        P_th_cgs,
        P_CR_cgs
        )

end


"""
    GadgetPhysical( DT::DataType, 
                    l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;
                    a_scale::Real=1.0, hpar::Real=1.0,
                    γ_th::Real=5.0/3.0, γ_CR::Real=4.0/3.0, xH::Real=0.76)

Set up a unit struct with a given `DataType`.
"""
GadgetPhysical( DT::DataType, 
                l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;
                a_scale::Real=1.0, hpar::Real=1.0,
                γ_th::Real=5.0/3.0, γ_CR::Real=4.0/3.0, xH::Real=0.76) = 
    GadgetPhysical( DT(l_unit), DT(m_unit), DT(v_unit),
                    a_scale=DT(a_scale), hpar=DT(hpar), 
                    γ_th=DT(γ_th), γ_CR=DT(γ_CR), xH=DT(xH) )