module AtomI

export Atoma

using StaticArrays

struct Atoma
    id::Int32
    atom_name::String
    alternativ_location::String
    parrent::Ref{}
    coords::SVector{3, Float32}
    occupancy::Float32
    bfactor::Float32
    segment_id::Int16
    element::String
    charge::Float32
    model::String
    chain::String
    compound::String
    idcompound::Int16
            
    properties::Dict{String, String}
    bonds::Vector{Ref}
end

get_ID(ob::Atoma) = ob.id
get_name(ob::Atoma) = ob.atom_name
get_alternativLoc(ob::Atoma) = ob.alternativ_location
get_parent(ob::Atoma) = parent
get_occupancy(ob::Atoma) = occupancy
get_location(ob::Atoma) = coords
get_BFactor(ob::Atoma) = bfactor
get_segmentID(ob::Atoma) = segment_id
get_element(ob::Atoma) = element
get_charge(ob::Atoma) = charge
get_property(ob::Atoma, propkey) = haskey(properties, propkey) ? properties[propkey] : missing
get_bonds(ob::Atoma) = bonds

set_parent(ob::Atoma, parent) = begin Atoma.parent = parent end
set_charge(ob::Atoma, charge) = begin Atoma.charge = charge end
set_property(ob::Atoma, property, value) = begin properties[propery] = value end
set_bond(ob::Atoma, new_bond) = begin push!(bonds, new_bond) end

end     # AtomI
