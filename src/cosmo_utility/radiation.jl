"""
    mJy_to_W( c::Cosmology.AbstractCosmology, S::Unitful.AbstractQuantity, z::Real )

Converts synchrotron flux density `S` in `[mJy]` to `[W/Hz]` for a given redshift `z` and cosmology `c`.
Conserves unit information.
"""
function mJy_to_W( c::Cosmology.AbstractCosmology, S::Unitful.AbstractQuantity, z::Real )

    # Luminsity Distance in Mpc
    dL = luminosity_dist(c, z)

    # P(ν)= S(ν) * 4π * dL^2
    # dL is the luminosity distance
    return S * 4π * dL^2 |> u"W/Hz"
end

"""
    mJy_to_W( c::Cosmology.AbstractCosmology, S::Real, z::Real )

Converts synchrotron flux density `S` in `[mJy]` to `[W/Hz]` for a given redshift `z` and cosmology `c`.
Deletes unit information.
"""
mJy_to_W( c::Cosmology.AbstractCosmology, S::Real, z::Real ) = mJy_to_W(c, S * 1.0u"mJy", z) |> ustrip

"""
    mJy_to_W(S::Union{Real, Unitful.AbstractQuantity}, h::AbstractGadgetHeader)

Converts synchrotron flux density `S` in `[mJy]` to `[W/Hz]` for a given redshift and cosmology taken from `AbstractGadgetHeader`.
"""
mJy_to_W(S::Union{Real, Unitful.AbstractQuantity}, h::AbstractGadgetHeader) = mJy_to_W( cosmology(h), S, h.z )

"""
    mJy_to_W(S::Union{Real,Unitful.AbstractQuantity}, d::Unitful.AbstractQuantity)

Converts synchrotron flux density `S` in `[mJy]` to `[W/Hz]` for a given distance `d`.
"""
mJy_to_W(S::Unitful.AbstractQuantity, d::Unitful.AbstractQuantity) = S * 4π * d^2 |> u"W/Hz"

"""
    mJy_to_W(S::Union{Real,Unitful.AbstractQuantity}, d::Real)

Converts synchrotron flux density `S` in `[mJy]` to `[W/Hz]` for a given distance `d` in `[Mpc]`.
"""
mJy_to_W(S::Real, d::Real) = mJy_to_W(S * 1.0u"mJy", d * 1.0u"Mpc") |> ustrip

"""
    mJy_to_W(S::Union{Real,Unitful.AbstractQuantity}, d::Real)

Converts synchrotron flux density `S` in `[mJy]` to `[W/Hz]` for a given distance `d` in `[Mpc]`.
"""
mJy_to_W(S::Unitful.AbstractQuantity, d::Real) = mJy_to_W(S, d * 1.0u"Mpc") |> ustrip


