module GadgetUnits

    export GadgetPhysical,
           GadgetPhysicalUnits,
           @u_str,
           strip_unit

    using Unitful
    using UnitfulAstro


    """
        GadgetPhysicalUnits(l_unit::Float64=3.085678e21, m_unit::Float64=1.989e43, v_unit::Float64=1.e5;
                            a_scale::Float64=1.0, hpar::Float64=1.0,
                            γ_th::Float64=5.0/3.0, γ_CR::Float64=4.0/3.0, xH::Float64=0.76)

    Creates a datatype GadgetPhysicalUnits which holds the conversion factors between comoving code units and physical units.
    Stores the unit information which can be converted with Unitful or UnitfulAstro.

    # Examples
    ```jldoctest
    julia> GU = GadgetPhysicalUnits()
    GadgetPhysicalUnits(3.085678e21 cm, 100000.0 cm s^-1, 1.989e43 g, 3.085678e16 s, 977.7923542981722 Myr, 1.989e53 erg, 1.2414361549102458e65 eV, 1.0 Gs, 6.769911178294544e-22 g cm^-3, 743179.9340255889 N_e/cm^3, 47.50882854026919 K, 6.769911178294544e-12 Ba, 6.769911178294544e-12 Ba)
    julia> GU.x_cgs
    3.085678e21 cm
    ```

    # Keyword Arguments
    - `a_scale::Float64 = 1.0`:  Cosmological scale factor of the simulation. Can be passed with the header `h` as `h.time`.
    - `hpar::Float64 = 1.0`:     Hubble constant as 'little h'. Can be passed with header `h` as `h.h0`.
    - `γ_th::Float64 = 5.0/3.0`: Adiabatic index of gas.
    - `γ_CR::Float64 = 4.0/3.0`: Adiabatic index of cosmic ray component.
    - `xH::Float64 = 0.76`:      Hydrogen fraction of the simulation, if run without chemical model.

    # Fields
    | Name                     | Meaning                        |
    |: ----------------------- |:-----------------------------  |
    | `x_cgs::Float64`         | position in cm                 |
    | `x_kpc::Float64`         | position in kpc                |
    | `v_cgs::Float64`         | velocity in cm/s               |
    | `v_kms::Float64`         | velocity in km/s               |
    | `m_cgs::Float64`         | mass in g                      |
    | `m_msun::Float64`        | mass in Msun                   |
    | `t_s::Float64`           | time in sec                    |
    | `t_Myr::Float64`         | time in Myr                    |
    | `E_cgs::Float64`         | energy in erg                  |
    | `E_eV::Float64`          | energy in eV                   |
    | `B_cgs::Float64`         | magnetic field in Gauss        |
    | `rho_cgs::Float64`       | density in ``g/cm^3``          |
    | `rho_ncm3::Float64`      | density in ``n_p/cm^3``        |
    | `T_K::Float64`           | temperature in K               |
    | `T_eV::Float64`          | temperature in eV              |
    | `P_th_cgs::Float64`      | thermal pressure in Ba         |
    | `P_CR_cgs::Float64`      | cosmic ray pressure in Ba      |

    """
    struct GadgetPhysicalUnits

        x_cgs::typeof(1.0u"cm")         # position in cm
        x_kpc::typeof(1.0u"kpc")        # position in kpc

        v_cgs::typeof(1.0u"cm/s")       # velocity in cm/s
        v_kms::typeof(1.0u"km/s")       # velocity in km/s

        m_cgs::typeof(1.0u"g")          # mass in g
        m_msun::typeof(1.0u"Msun")      # mass in Msun

        t_s::typeof(1.0u"s")            # time in sec
        t_Myr::typeof(1.0u"Myr")        # time in Myr

        E_cgs::typeof(1.0u"erg")        # energy in erg
        E_eV::typeof(1.0u"eV")          # energy in eV

        B_cgs::typeof(1.0u"Gs")         # magnetic field in Gauss

        rho_cgs::typeof(1.0u"g/cm^3")   # density in g/cm^3
        rho_ncm3::typeof(1.0u"cm^-3")     # density in mp/cm^3

        T_K::typeof(1.0u"K")            # temperature in K
        T_eV::typeof(1.0u"eV")            # temperature in eV

        P_th_cgs::typeof(1.0u"erg/cm^3")      # thermal pressure in Ba
        P_CR_cgs::typeof(1.0u"erg/cm^3")      # cosmic ray pressure in Ba

        function GadgetPhysicalUnits(l_unit::Float64=3.085678e21, m_unit::Float64=1.989e43, v_unit::Float64=1.e5;
                                    a_scale::Float64=1.0, hpar::Float64=1.0,
                                    γ_th::Float64=5.0/3.0, γ_CR::Float64=4.0/3.0, xH::Float64=0.76)

            # Gadget units are given in cgs
            l_unit *= 1.0u"cm"
            m_unit *= 1.0u"g"
            v_unit *= 1.0u"cm/s"

            # convert comoving output to physical units
            x_cgs   = l_unit * a_scale / hpar
            x_kpc   = x_cgs |> u"kpc"

            v_cgs   = v_unit * sqrt(a_scale)
            v_kms   = v_cgs |> u"km/s"

            m_cgs   = m_unit / hpar
            m_msun  = m_cgs |> u"Msun"

            t_unit  = l_unit / v_unit
            t_s     = t_unit * sqrt(a_scale) / hpar  # in sec
            t_Myr   = t_s |> u"Myr"

            E_cgs = m_cgs * v_cgs^2 |> u"erg"
            E_eV = E_cgs |> u"eV"

            B_cgs = 1.0u"Gs"    # gadget outputs in cgs

            rho_cgs = m_unit/l_unit^3 * hpar^2 / a_scale^3

            yhelium = ( 1.0 - xH ) / ( 4.0 * xH )
            mean_mol_weight = (1.0 + 4.0 * yhelium) / (1.0 + 3.0 * yhelium + 1.0)

            n2ne =  (  xH + 0.5 * ( 1.0 - xH ) ) / 
                    ( 2.0 * xH + 0.75 * ( 1.0 - xH ) )       # conversion n_pat -> n_electrons
            umu  =  4.0 / ( 5.0 * xH + 3.0 )                 # mean molucular weight in hydr. mass

            rho_ncm3 = rho_cgs * n2ne/( umu * 1.0u"mp" ) |> u"cm^-3"
            #rho_ncm3 = rho_cgs |> u"n_p"

            yhelium = ( 1.0 - xH ) / ( 4.0 * xH )
            mean_mol_weight = (1.0 + 4.0 * yhelium) / (1.0 + 3.0 * yhelium + 1.0)

            T_cgs = (γ_th - 1.0) * v_unit^2 * 1.0u"mp" * mean_mol_weight / 1.0u"k" |> u"K"
            T_eV  = T_cgs * 1.0u"k" |> u"eV"

            P_th_cgs = a_scale^(-3) * E_cgs / l_unit^3 * hpar^2  |> u"erg/cm^3"
            P_CR_cgs = a_scale^(-4) * E_cgs / l_unit^3 * hpar^2  |> u"erg/cm^3"

            new(x_cgs, x_kpc,
                v_cgs, v_kms,
                m_cgs, m_msun,
                t_s, t_Myr,
                E_cgs, E_eV,
                B_cgs,
                rho_cgs, rho_ncm3,
                T_cgs, T_eV,
                P_th_cgs,
                P_CR_cgs
                )

        end
    end


    """
        GadgetPhysical(l_unit::Float64=3.085678e21, m_unit::Float64=1.989e43, v_unit::Float64=1.e5;
                    a_scale::Float64=1.0, hpar::Float64=1.0,
                    γ_th::Float64=5.0/3.0, γ_CR::Float64=4.0/3.0, xH::Float64=0.76)

    Creates a datatype GadgetPhysical which holds the conversion factors between comoving code units and physical units, without unit information.

    # Examples
    ```jldoctest
    julia> GU = GadgetPhysical()
    GadgetPhysicalUnits(3.085678e21, 100000.0, 1.989e43, 3.085678e16, 977.7923542981722, 1.989e53, 1.2414361549102458e65, 1.0, 6.769911178294544e-22, 743179.9340255889, 47.50882854026919, 6.769911178294544e-12, 6.769911178294544e-12)
    julia> GU.x_cgs
    3.085678e21
    ```

    Keyword arugments specify:
    # Arguments
    - `a_scale::Float64 = 1.0`:  Cosmological scale factor of the simulation. Can be passed with the header `h` as `h.time`.
    - `hpar::Float64 = 1.0`:     Hubble constant as 'little h'. Can be passed with header `h` as `h.h0`.
    - `γ_th::Float64 = 5.0/3.0`: Adiabatic index of gas.
    - `γ_CR::Float64 = 4.0/3.0`: Adiabatic index of cosmic ray component.
    - `xH::Float64 = 0.76`:      Hydrogen fraction of the simulation, if run without chemical model.

    # Fields
    | Name                     | Meaning                        |
    |: ----------------------- |:-----------------------------  |
    | `x_cgs::Float64`         | position in cm                 |
    | `x_kpc::Float64`         | position in kpc                |
    | `v_cgs::Float64`         | velocity in cm/s               |
    | `v_kms::Float64`         | velocity in km/s               |
    | `m_cgs::Float64`         | mass in g                      |
    | `m_msun::Float64`        | mass in Msun                   |
    | `t_s::Float64`           | time in sec                    |
    | `t_Myr::Float64`         | time in Myr                    |
    | `E_cgs::Float64`         | energy in erg                  |
    | `E_eV::Float64`          | energy in eV                   |
    | `B_cgs::Float64`         | magnetic field in Gauss        |
    | `rho_cgs::Float64`       | density in ``g/cm^3``          |
    | `rho_ncm3::Float64`      | density in ``n_p/cm^3``        |
    | `T_K::Float64`           | temperature in K               |
    | `T_eV::Float64`          | temperature in eV              |
    | `P_th_cgs::Float64`      | thermal pressure in Ba         |
    | `P_CR_cgs::Float64`      | cosmic ray pressure in Ba      |

    """
    struct GadgetPhysical

        x_cgs::Float64         # position in cm
        x_kpc::Float64         # position in kpc

        v_cgs::Float64         # velocity in cm/s
        v_kms::Float64         # velocity in km/s

        m_cgs::Float64         # mass in g
        m_msun::Float64        # mass in Msun

        t_s::Float64           # time in sec
        t_Myr::Float64         # time in Myr

        E_cgs::Float64         # energy in erg
        E_eV::Float64          # energy in eV

        B_cgs::Float64         # magnetic field in Gauss

        rho_cgs::Float64       # density in g/cm^3
        rho_ncm3::Float64      # density in N_p/cm^3

        T_K::Float64           # temperature in K
        T_eV::Float64          # temperature in eV

        P_th_cgs::Float64      # thermal pressure in Ba
        P_CR_cgs::Float64      # cosmic ray pressure in Ba

        function GadgetPhysical(l_unit::Float64=3.085678e21, m_unit::Float64=1.989e43, v_unit::Float64=1.e5;
                                a_scale::Float64=1.0, hpar::Float64=1.0,
                                γ_th::Float64=5.0/3.0, γ_CR::Float64=4.0/3.0, xH::Float64=0.76)

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
            x_cgs   = l_unit * a_scale / hpar
            v_cgs   = v_unit * sqrt(a_scale)
            m_cgs   = m_unit / hpar
            t_unit  = l_unit / v_unit
            t_s     = t_unit * sqrt(a_scale) / hpar  # in sec
            t_Myr   = t_s / 3.15576e13

            E_cgs = m_cgs * v_cgs^2
            E_eV = E_cgs

            B_cgs = 1.0   # gadget outputs in cgs

            rho_cgs = m_unit/l_unit^3 * hpar^2 / a_scale^3
            rho_ncm3 = rho_cgs * n2ne/( umu * mp )

        
            T_cgs = (γ_th - 1.0) * v_unit^2 * mean_mol_weight * mp / kB
            T_eV  = T_cgs * kB / eV2cgs

            P_th_cgs = a_scale^(-3) * E_cgs / l_unit^3 * hpar^2
            P_CR_cgs = a_scale^(-4) * E_cgs / l_unit^3 * hpar^2

            new(x_cgs, x_cgs/3.085678e21, 
                v_cgs, v_cgs*1.e-5,
                m_cgs, m_cgs/1.989e43,
                t_s, t_Myr,
                E_cgs, E_eV,
                B_cgs,
                rho_cgs, rho_ncm3,
                T_cgs, T_eV,
                P_th_cgs,
                P_CR_cgs
                )

        end
    end

    """
        strip_unit(a)

    Strips unit information from the input to convert back to simple `AbstractFloat` datatypes.

    # Examples
    ```jldoctest
    julia> x = 3.085678e21u"cm"
    3.085678e21 cm
    julia> strip_unit(x)
    3.085678e21
    ```

    """
    @inline function strip_unit(a)
        return a / unit(a)
    end

end # module
