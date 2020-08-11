using Documenter
using GadgetUnits

makedocs(
    sitename="GadgetUnits",
    format=Documenter.HTML(),
    modules=[GadgetUnits]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/LudwigBoess/GadgetUnits.jl.git"
)