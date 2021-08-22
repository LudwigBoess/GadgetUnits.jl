import Cosmology.cosmology 

"""
    cosmology(h::SnapshotHeader)

Defines a `Cosmology.AbstractCosmology` from the properties of `h`.
"""
cosmology(h::SnapshotHeader) = cosmology(h=h.h0, OmegaM=h.omega_0)