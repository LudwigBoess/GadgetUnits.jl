import Cosmology.age
import Cosmology.lookback_time
using Roots

"""
    age(h::AbstractGadgetHeader, units::Bool=false)

Computes the age of the universe given the properties from `AbstractGadgetHeader`.
"""
function age(h::AbstractGadgetHeader, units::Bool=false)
    
    c = cosmology(h)
    t = age(c, h.z)

    if units
        return t
    else
        return t |> ustrip
    end
end

"""
    lookback_time(h::AbstractGadgetHeader, units::Bool=false)

Computes the lookback tima given the properties from `AbstractGadgetHeader`.
"""
function lookback_time(h::AbstractGadgetHeader, units::Bool=false)

    c = cosmology(h)
    t = lookback_time(c, h.z)

    if units
        return t
    else
        return t |> ustrip
    end
end



"""
    redshift(c::AbstractCosmology, t_Gyrs::Real)

Computes the redshift at a given age of the universe in Gyrs for a cosmology `t`.
"""
function redshift(t_Gyrs::Real, c::Cosmology.AbstractCosmology)

    age_helper(z) = ustrip(age(c, z)) - t_Gyrs

    return find_zero(age_helper, 0.0 )
end



"""
    redshift(t::Unitful.AbstractQuantity, h::AbstractGadgetHeader)

See [`redshift`](@ref).
"""
redshift(t::Unitful.AbstractQuantity, h::AbstractGadgetHeader) = redshift( ustrip( t |> u"Gyr" ), cosmology(h))


"""
    redshift(t::Real, h::AbstractGadgetHeader)

See [`redshift`](@ref).
"""
redshift(t::Real, h::AbstractGadgetHeader) = redshift( t, cosmology(h))


"""
    redshift(t::Unitful.AbstractQuantity, h::AbstractGadgetHeader)

See [`redshift`](@ref).
"""
redshift(t::Unitful.AbstractQuantity, c::Cosmology.AbstractCosmology) = redshift( ustrip( t |> u"Gyr" ), c )