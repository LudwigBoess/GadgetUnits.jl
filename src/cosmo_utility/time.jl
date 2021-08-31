import Cosmology.age 

"""
    age(h::AbstractGadgetHeader)

Computes the age of the universe given the properties from `AbstractGadgetHeader`.
"""
function age(h::AbstractGadgetHeader, units::Bool=true)
    
    c = cosmology(h)
    t = age(c, h.z)

    if units
        return t
    else
        return t |> ustrip
    end
end