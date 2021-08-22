import Cosmology.age 

"""
    age(h::SnapshotHeader)

Computes the age of the universe given the properties from `SnapshotHeader`.
"""
function age(h::SnapshotHeader, units::Bool=true)
    
    c = cosmology(h)
    t = age(c, h.z)

    if units
        return t
    else
        return t |> ustrip
    end
end