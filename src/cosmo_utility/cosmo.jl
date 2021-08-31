import Cosmology.cosmology 

"""
    cosmology(h::AbstractGadgetHeader)

Defines a `Cosmology.AbstractCosmology` from the properties of `h`.
"""
cosmology(h::AbstractGadgetHeader) = cosmology(h=h.h0, OmegaM=h.omega_0)