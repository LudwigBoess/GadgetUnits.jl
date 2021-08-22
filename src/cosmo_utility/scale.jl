"""
    arcmin_to_kpc(c::Cosmology.AbstractCosmology, θ::Unitful.AbstractQuantity, z::Real)

Convert `arcmin` to `kpc` for a given redshift `z` and cosmology `c`.
Conserves unit information.
"""
function arcmin_to_kpc(c::Cosmology.AbstractCosmology, θ::Unitful.AbstractQuantity, z::Real)

    # angular diameter distance in Mpc
    dA = angular_diameter_dist(c, z)

    return dA * 2π / 360 / 60.0u"arcminute" * θ |> u"kpc"
end

"""
    arcmin_to_kpc(c::Cosmology.AbstractCosmology, θ::Real, z::Real)

Convert `arcmin` to `kpc` for a given redshift `z` and cosmology `c`.
Deletes unit information.
"""
arcmin_to_kpc(c::Cosmology.AbstractCosmology, θ::Real, z::Real ) = arcmin_to_kpc( cosmology(h=h.h0, OmegaM=h.omega_0), θ * 1.0u"arcminute", h.z) |> ustrip

"""
    arcmin_to_kpc(θ::Real, h::SnapshotHeader)

Convert `arcmin` to `kpc` for a given redshift and cosmology taken from `SnapshotHeader`.
"""
arcmin_to_kpc(θ::Real, h::SnapshotHeader) = arcmin_to_kpc( cosmology(h=h.h0, OmegaM=h.omega_0), θ, h.z) |> ustrip

"""
    arcmin_to_kpc(θ::Unitful.AbstractQuantity, h::SnapshotHeader)

Convert `arcmin` to `kpc` for a given redshift and cosmology taken from `SnapshotHeader`.
"""
arcmin_to_kpc(θ::Unitful.AbstractQuantity, h::SnapshotHeader) = arcmin_to_kpc( cosmology(h=h.h0, OmegaM=h.omega_0), θ, h.z)