module GadgetUnits

    export GadgetPhysical,
           GadgetPhysicalUnits,
           strip_unit

    using Unitful
    using UnitfulAstro


    """
        GadgetPhysicalUnits(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;
                            a_scale::T=1.0, hpar::T=1.0,
                            γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T

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
    - `a_scale::T = 1.0`:  Cosmological scale factor of the simulation. Can be passed with the header `h` as `h.time`.
    - `hpar::T = 1.0`:     Hubble constant as 'little h'. Can be passed with header `h` as `h.h0`.
    - `γ_th::T = 5.0/3.0`: Adiabatic index of gas.
    - `γ_CR::T = 4.0/3.0`: Adiabatic index of cosmic ray component.
    - `xH::T = 0.76`:      Hydrogen fraction of the simulation, if run without chemical model.

    # Fields
    | Name                         | Meaning                        |
    |: --------------------------- |:-----------------------------  |
    | `x_cgs::Quantity{T}`         | position in cm                 |
    | `x_kpc::Quantity{T}`         | position in kpc                |
    | `v_cgs::Quantity{T}`         | velocity in cm/s               |
    | `v_kms::Quantity{T}`         | velocity in km/s               |
    | `m_cgs::Quantity{T}`         | mass in g                      |
    | `m_msun::Quantity{T}`        | mass in Msun                   |
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

    """
    struct GadgetPhysicalUnits{T}

        x_cgs::Quantity{T}        # position in cm
        x_kpc::Quantity{T}        # position in kpc

        v_cgs::Quantity{T}     # velocity in cm/s
        v_kms::Quantity{T}     # velocity in km/s

        m_cgs::Quantity{T}     # mass in g
        m_msun::Quantity{T}    # mass in Msun

        t_s::Quantity{T}          # time in sec
        t_Myr::Quantity{T}        # time in Myr

        E_cgs::Quantity{T}        # energy in erg
        E_eV::Quantity{T}         # energy in eV

        B_cgs::Quantity{T}     # magnetic field in Gauss

        rho_cgs::Quantity{T}   # density in g/cm^3
        rho_ncm3::Quantity{T}    # density in mp/cm^3

        T_K::Quantity{T}           # temperature in K
        T_eV::Quantity{T}            # temperature in eV

        P_th_cgs::Quantity{T}      # thermal pressure in Ba
        P_CR_cgs::Quantity{T}      # cosmic ray pressure in Ba

        function GadgetPhysicalUnits(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;
                                    a_scale::T=1.0, hpar::T=1.0,
                                    γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T

            # Gadget units are given in cgs
            l_unit *= one(T)*u"cm"
            m_unit *= one(T)*u"g"
            v_unit *= one(T)*u"cm/s"

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

            new{T}(x_cgs, x_kpc,
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
        GadgetPhysicalUnits( DT::DataType, 
                             l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;
                             a_scale::Real=1.0, hpar::Real=1.0,
                             γ_th::Real=5.0/3.0, γ_CR::Real=4.0/3.0, xH::Real=0.76)

    Set up a unit struct with a given `DataType`.
    """
    GadgetPhysicalUnits( DT::DataType, 
                         l_unit::Real=3.085678e21, m_unit::Real=1.989e43, v_unit::Real=1.e5;
                         a_scale::Real=1.0, hpar::Real=1.0,
                         γ_th::Real=5.0/3.0, γ_CR::Real=4.0/3.0, xH::Real=0.76) = 
        GadgetPhysicalUnits(DT(l_unit), DT(m_unit), DT(v_unit),
                            a_scale=DT(a_scale), hpar=DT(hpar), 
                            γ_th=DT(γ_th), γ_CR=DT(γ_CR), xH=DT(xH))


    """
        GadgetPhysical(l_unit::T=3.085678e21, m_unit::T=1.989e43, v_unit::T=1.e5;
                    a_scale::T=1.0, hpar::T=1.0,
                    γ_th::T=5.0/3.0, γ_CR::T=4.0/3.0, xH::T=0.76) where T

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
    - `a_scale::T = 1.0`:  Cosmological scale factor of the simulation. Can be passed with the header `h` as `h.time`.
    - `hpar::T = 1.0`:     Hubble constant as 'little h'. Can be passed with header `h` as `h.h0`.
    - `γ_th::T = 5.0/3.0`: Adiabatic index of gas.
    - `γ_CR::T = 4.0/3.0`: Adiabatic index of cosmic ray component.
    - `xH::T = 0.76`:      Hydrogen fraction of the simulation, if run without chemical model.

    # Fields
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

    """
    struct GadgetPhysical{T}

        x_cgs::T         # position in cm
        x_kpc::T         # position in kpc

        v_cgs::T         # velocity in cm/s
        v_kms::T         # velocity in km/s

        m_cgs::T         # mass in g
        m_msun::T        # mass in Msun

        t_s::T           # time in sec
        t_Myr::T         # time in Myr

        E_cgs::T         # energy in erg
        E_eV::T          # energy in eV

        B_cgs::T         # magnetic field in Gauss

        rho_cgs::T       # density in g/cm^3
        rho_ncm3::T      # density in N_p/cm^3

        T_K::T           # temperature in K
        T_eV::T          # temperature in eV

        P_th_cgs::T      # thermal pressure in Ba
        P_CR_cgs::T      # cosmic ray pressure in Ba

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

            new{T}(x_cgs, x_cgs/T(3.085678e21), 
                v_cgs, v_cgs*T(1.e-5),
                m_cgs, m_cgs/T(1.989e33),
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
