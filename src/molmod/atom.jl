# module AtomI

export Atoma
export StructModel
export AtomsGroup
export PDBsChain

using StaticArrays


#= struct Atoma
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
 =#

struct Atoma
    group_PDB::Symbol                               #1  HETATOM, ATOM
    id::Int32                                       #2  Atom index
    type_symbol::Symbol                             #3  Element
    label_atom_id::Symbol                           #4  Type Atom in compound(auth_atom_id)
    label_alt_id::Symbol                            #5  Alternative in Structure
    label_comp_id::Symbol                           #6  Compound(atomic group) (auth_comp_id)
    label_asym_id::Symbol                           #7  Chain Id (auth_asym_id)
    label_entity_id::Symbol                         #8  ?
    label_seq_id::Int32                             #9  Compound(atomic group) index (auth_seq_id)
    pdbx_PDB_ins_code::Symbol                       #10 ?
    Cartn_x                                         #11 Coords
    Cartn_y                                         #12
    Cartn_z                                         #13
    occupancy::Float32                              #14 Part position detection
    B_iso_or_equiv::Float32                         #15 B-factor
    pdbx_formal_charge::Union{Float32, Nothing}     #16 Formal pdb charge
    auth_seq_id::Int32                              #17 Compound(atomic group) index (label_seq_id)
    auth_comp_id::Symbol                            #18 Compound(atomic group)(label_comp_id)
    auth_asym_id::Symbol                            #19 Chain Id (label_asym_id)
    auth_atom_id::Symbol                            #20 Type Atom in compound(label_atom_id)
    pdbx_PDB_model_num::Int32                       #21 Model Index
    parents::Dict{Type, Ref{}}                      # links to struct hierarhy elements
    function Atoma(ar::Vector{SubString{String}})
        #println(ar)
        if length(ar)!=21 error("atomic notation of incorrect length") end
        latom_id = (ar[4]  in ("?", ".") ? Symbol(ar[20]) : Symbol(ar[4]))
        aatom_id = (ar[20] in ("?", ".") ? Symbol(ar[4]) : Symbol(ar[20]))
        lcomp_id = (ar[6]  in ("?", ".") ? Symbol(ar[18]) : Symbol(ar[6]))
        acomp_id = (ar[18] in ("?", ".") ? Symbol(ar[6]) : Symbol(ar[18]))
        lasym_id = (ar[7]  in ("?", ".") ? Symbol(ar[19]) : Symbol(ar[7]))
        aasym_id = (ar[19] in ("?", ".") ? Symbol(ar[7]) : Symbol(ar[19]))
        lseq_id  = (ar[9]  in ("?", ".") ? parse(Int32, ar[17]) : parse(Int32, ar[9]))
        aseq_id  = (ar[17] in ("?", ".") ? parse(Int32, ar[9]) : parse(Int32, ar[17]))
        new(Symbol(ar[1]), 
            parse(Int32, ar[2]), 
            Symbol(ar[3]), 
            latom_id, 
            Symbol(ar[5]),
            lcomp_id,
            lasym_id,
            Symbol(ar[8]),
            lseq_id,
            Symbol(ar[10]),
            parse(Float32, ar[11]),
            parse(Float32, ar[12]),
            parse(Float32, ar[13]),
            parse(Float32, ar[14]),
            parse(Float32, ar[15]),
            tryparse(Float32, ar[16]),
            aseq_id,
            acomp_id,
            aasym_id,
            aatom_id,
            parse(Int32, ar[21]),
            Dict{Type, Ref{}}())
    end
end

struct AtomsGroup
    id::Tuple{Int32, Symbol, Int32}
    compoundname::Symbol
    childs::Dict{Type, Vector{Ref{}}}
    parent::Dict{Type, Ref{}}
    #propertyes::Dict{AbstractString, AbstractString}
end

struct PDBsChain
    id::Tuple{Int32, Symbol}
    childs::Dict{Type, Vector{Ref{}}}
    parent::Dict{Type, Ref{}}
end

struct StructModel
    id::Int32
    childs::Dict{Type, Vector{Ref{}}}
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
