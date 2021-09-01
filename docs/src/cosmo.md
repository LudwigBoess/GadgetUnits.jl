```@meta
CurrentModule = GadgetUnits
DocTestSetup = quote
    using GadgetUnits
end
```

# Cosmological conversions

`GadgetUnits.jl` supplied some conversion functions to convert observables in given cosmologies. These functions rely on [Cosmology.jl](https://github.com/JuliaAstro/Cosmology.jl), specifically a struct of `AbstractCosmology`.

## Cosmology

To define an `AbstractCosmology` from the properties of the simulation you can use the multiple dispatch function
```julia
c = cosmology(h::AbstractGadgetHeader)
```

## Time

To get the current age of the universe you can use [`age`](@ref)

```julia
t = age(h::AbstractGadgetHeader, units::Bool=true)
```

If the optional parameter `units` is set to `true` it will preserve the unit property. Otherwise it will return a `Real` in units `Gyrs`.

## Scale

You can convert arcminutes to kpc for an object at redshift `z` with [`armin_to_kpc`](@ref).
This can be used either with a pre-defined `Cosmology.AbstractCosmology` 
```julia
arcmin_to_kpc(c::Cosmology.AbstractCosmology, θ::Real, z::Real )
```
or can use a `AbstractGadgetHeader` struct from [`GadgetIO.jl`](https://github.com/LudwigBoess/GadgetIO.jl) which defines the `AbstractCosmology` on the fly:

```julia
arcmin_to_kpc(θ::Real, h::AbstractGadgetHeader)
```

If you input `θ::Real` the function will return `kpc`, but remove the unit property.
Optionally you can define `θ::Unitful.AbstractQuantity` which preserves the unit property.

## Radiation

To convert synchrotron radiation from `mJy` to `W/Hz` you can use [`mJy_to_W`](@ref).
Like in the case of [`Scale`](@ref) you can use either a pre-defined `AbstractCosmology`
```julia
mJy_to_W( c::Cosmology.AbstractCosmology, S::Real, z::Real )
```

or use a `AbstractGadgetHeader`

```julia
mJy_to_W(S::Union{Real, Unitful.AbstractQuantity}, h::AbstractGadgetHeader)
```