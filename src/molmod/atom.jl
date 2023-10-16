# module AtomI

export Atoma
export StructModel
export AtomsGroup
export PDBsChain

using StaticArrays

struct Atoma
    id::Union{Int32, Nothing}
    hetatom::Bool
    atom_name::Union{String, Nothing}
    alternativ::Union{String, Nothing}
    coords::SVector{3, T} where T<:Union{Float32, Nothing}
    occupancy::Union{Float32, Nothing}
    bfactor::Union{Float32, Nothing}

    element::Union{String, Nothing}
    charge::Union{Float32, Nothing}
    
    model::Union{Int32, Nothing}
    chain::Union{String, Nothing}
    compound::Union{String, Nothing}
    idcompound::Union{Int32, Nothing}

    i_categori::Vector{Int8}
    parents::Dict{String, Ref{}}
    #properties::Dict{String, Union{Number, AbstractString}}
    #bonds::Vector{Ref}
end

struct AtomsGroup
    id::Tuple{Int32, AbstractString, Int32}
    compoundname::String
    childs::Dict{String, Vector{Ref{}}}
    parent::Dict{String, Ref{}}
    #propertyes::Dict{AbstractString, AbstractString}
end

struct PDBsChain
    id::Tuple{Int32, AbstractString}
    childs::Dict{String, Vector{Ref{}}}
    parent::Dict{String, Ref{}}
end

struct StructModel
    id::Int32
    childs::Dict{String, Vector{Ref{}}}
end


# get_ID(ob::Atoma) = ob.id
# get_name(ob::Atoma) = ob.atom_name
# get_alternativLoc(ob::Atoma) = ob.alternativ_location
# get_parent(ob::Atoma) = parent
# get_occupancy(ob::Atoma) = occupancy
# get_location(ob::Atoma) = coords
# get_BFactor(ob::Atoma) = bfactor
# get_segmentID(ob::Atoma) = segment_id
# get_element(ob::Atoma) = element
# get_charge(ob::Atoma) = charge
# get_property(ob::Atoma, propkey) = haskey(properties, propkey) ? properties[propkey] : missing
# get_bonds(ob::Atoma) = bonds

# set_parent(ob::Atoma, parent) = begin Atoma.parent = parent end
# set_charge(ob::Atoma, charge) = begin Atoma.charge = charge end
# set_property(ob::Atoma, property, value) = begin properties[propery] = value end
# set_bond(ob::Atoma, new_bond) = begin push!(bonds, new_bond) end

# end     # AtomI
