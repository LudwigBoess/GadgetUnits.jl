```@meta
CurrentModule = GadgetUnits
DocTestSetup = quote
    using GadgetUnits
end
```

# GadgetUnits.jl

This package provides some basic unit conversion functionality to work with the SPH code [Gadget2](https://wwwmpa.mpa-garching.mpg.de/gadget/) by Volker Springel.
It is based on the unit conversions listed by [Klaus Dolag](https://www.usm.uni-muenchen.de/~dolag/) on his HowTo page (restricted).
Further functionality provides tools for cosmological simulations.

If you use `GadgetUnits.jl` in publications please [cite](https://zenodo.org/badge/latestdoi/277262050) it.

# Table of Contents

```@contents
Pages = [ "index.md",
          "install.md",
          "conversion_structs.md",
          "cosmo.md",
          "api.md" 
        ]
Depth = 3
```