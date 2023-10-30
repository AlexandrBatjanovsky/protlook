# module AtomI

export Atoma
export StructModel
export AtomsGroup
export PDBsChain

using StaticArrays: SVector

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
    Cartn_x::Float32                                #11 Coords
    Cartn_y::Float32                                #12
    Cartn_z::Float32                                #13
    occupancy::Float32                              #14 Part position detection
    B_iso_or_equiv::Float32                         #15 B-factor
    pdbx_formal_charge::Union{Float32, Nothing}     #16 Formal pdb charge
    auth_seq_id::Int32                              #17 Compound(atomic group) index (label_seq_id)
    auth_comp_id::Symbol                            #18 Compound(atomic group)(label_comp_id)
    auth_asym_id::Symbol                            #19 Chain Id (label_asym_id)
    auth_atom_id::Symbol                            #20 Type Atom in compound(label_atom_id)
    pdbx_PDB_model_num::Int32                       #21 Model Index
    calc_flag::Union{Int16, Nothing}                 # calc or intense 
    Cartn_x_esd::Union{Float32, Nothing}
    Cartn_y_esd::Union{Float32, Nothing}
    Cartn_z_esd::Union{Float32, Nothing}
    occupancy_esd::Union{Float32, Nothing}
    B_iso_or_equiv_esd::Union{Float32, Nothing}
    parents::Dict{Type, Ref{}}                      # links to struct hierarhy elements
    function Atoma(ar::NamedTuple)
        latom_id = ar.label_atom_id  in ("?", ".") ? Symbol(ar.auth_atom_id) : Symbol(ar.label_atom_id)
        aatom_id = ar.auth_atom_id in ("?", ".") ? Symbol(ar.label_atom_id) : Symbol(ar.auth_atom_id)
        lcomp_id = ar.label_comp_id  in ("?", ".") ? Symbol(ar.auth_comp_id) : Symbol(ar.label_comp_id)
        acomp_id = ar.auth_comp_id in ("?", ".") ? Symbol(ar.label_comp_id) : Symbol(ar.auth_comp_id)
        lasym_id = ar.label_asym_id  in ("?", ".") ? Symbol(ar.auth_asym_id) : Symbol(ar.label_asym_id)
        aasym_id = ar.auth_asym_id in ("?", ".") ? Symbol(ar.label_asym_id) : Symbol(ar.auth_asym_id)
        lseq_id  = ar.label_seq_id  in ("?", ".") ? parse(Int32, ar.auth_seq_id) : parse(Int32, ar.label_seq_id)
        aseq_id  = ar.auth_seq_id in ("?", ".") ? parse(Int32, ar.label_seq_id) : parse(Int32, ar.auth_seq_id)
        calc_flag = :calc_flag in keys(ar) ? tryparse(Int16, ar.calc_flag) : nothing
        Cartn_x_esd = :Cartn_x_esd in keys(ar) ? tryparse(Float32, ar.Cartn_x_esd) : nothing
        Cartn_y_esd = :Cartn_y_esd in keys(ar) ? tryparse(Float32, ar.Cartn_y_esd) : nothing
        Cartn_z_esd = :Cartn_z_esd in keys(ar) ? tryparse(Float32, ar.Cartn_z_esd) : nothing
        occupancy_esd = :occupancy_esd in keys(ar) ? tryparse(Float32, ar.occupancy_esd) : nothing
        B_iso_or_equiv_esd = :B_iso_or_equiv_esd in keys(ar) ? tryparse(Float32, ar.B_iso_or_equiv_esd) : nothing
        new(Symbol(ar.group_PDB), 
            parse(Int32, ar.id), 
            Symbol(ar.type_symbol), 
            latom_id, 
            Symbol(ar.label_alt_id),
            lcomp_id,
            lasym_id,
            Symbol(ar.label_entity_id),
            lseq_id,
            Symbol(ar.pdbx_PDB_ins_code),
            parse(Float32, ar.Cartn_x),
            parse(Float32, ar.Cartn_y),
            parse(Float32, ar.Cartn_z),
            parse(Float32, ar.occupancy),
            parse(Float32, ar.B_iso_or_equiv),
            tryparse(Float32, ar.pdbx_formal_charge),
            aseq_id,
            acomp_id,
            aasym_id,
            aatom_id,
            parse(Int32, ar.pdbx_PDB_model_num),
            calc_flag,
            Cartn_x_esd,
            Cartn_y_esd,
            Cartn_z_esd,
            occupancy_esd,
            B_iso_or_equiv_esd,
            Dict{Type, Ref{}}())
    end
end

struct Bondc
    comp_id::Symbol
    atom_id_1::Symbol
    atom_id_2::Symbol
    value_order::Symbol
    pdbx_aromatic_flag::Bool
    pdbx_stereo_config::Bool
    pdbx_ordinal::Int16
    function Bondc(ar::Vector{SubString{String}})
        new(Symbol(ar[1]),
            Symbol(ar[2]),
            Symbol(ar[3]),
            Symbol(ar[4]),
            ar[5] == "Y",
            ar[6] == "Y",
            parse(Int16, ar[7]))
    end
end

struct Atomc
    comp_id::Symbol
    atom_id::Symbol
    alt_atom_id::Symbol
    type_symbol::Symbol
    charge::Float32
    pdbx_align::Int16
    pdbx_aromatic_flag::Bool
    pdbx_leaving_atom_flag::Bool
    pdbx_stereo_config::Bool
    model_Cartn_x::Float32 
    model_Cartn_y::Float32
    model_Cartn_z::Float32
    pdbx_model_Cartn_x_ideal::Float32
    pdbx_model_Cartn_y_ideal::Float32
    pdbx_model_Cartn_z_ideal::Float32
    pdbx_component_atom_id::Symbol
    pdbx_component_comp_id::Symbol
    pdbx_ordinal::Int16
    # bonds::Dict{Symbol, Tuple{Bondc, Float32, Dict{Symbol, Float32}}}
    function Atomc(ar::Vector{SubString{String}})
        new(Symbol(ar[1]),
            Symbol(ar[2]),
            Symbol(ar[3]),
            Symbol(ar[4]),
            parse(Float32, ar[5]),
            parse(Int16, ar[6]),
            ar[7] == "Y",
            ar[8] == "Y",
            ar[9] == "Y",
            parse(Float32, ar[10]),
            parse(Float32, ar[11]),
            parse(Float32, ar[12]),
            parse(Float32, ar[13]),
            parse(Float32, ar[14]),
            parse(Float32, ar[15]),
            Symbol(ar[16]),
            Symbol(ar[17]),
            parse(Int16, ar[18]))
            # Dict{Symbol, Tuple{Bondc, Float32, Dict{Symbol, Float32}}}())
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
