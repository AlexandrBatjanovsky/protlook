module StereoChem


export loadCompounds, CmpAtomic, CompositionConfirmation

using ..Datas: settings
using ..Datas: Atoma, StructModel, AtomsGroup, PDBsChain
using ..Datas: Atomc, Bondc

using Distances

using Logging
using JLD2
using StaticArrays: SVector

SetCompositions = Dict{Symbol, @NamedTuple{atomc::Dict{Symbol, Atomc},
                                           catec::Dict{Symbol, Dict{Symbol, Vector{String}}}}}()

"""
    loadCompounds(CompoundNames)

Selects from components.cif a set of elementary compounds included in the database under study (amino acid residues, heterocompounds)

CompoundNames - Set of Symbols

#Examples
```julia-repl
julia> loadCompounds!(Set([:ALA, :ARG, :PRO]))
```
"""
function loadCompounds(compoundnames::Set{Symbol})
    global SetCompositions
    compoundlibfile = open(settings[:CmpLib], "r")
    compoundoffset = load(settings[:CmpInd], "compoundoffset")
    while !isempty(compoundnames)
        compoundname = pop!(compoundnames)
        # using indexes of position for compound in file
        seek(compoundlibfile, compoundoffset[compoundname])
        rl = readline(compoundlibfile)
        categories = Dict{Symbol, Dict{Symbol, Vector{String}}}()
        global curcategory = :absence
        atomic = Dict{Symbol, Atomc}()
        # two set of coordinates in components.cif with ? in all
        coordcm = Vector{SVector{3, Float32}}(); coordci = Vector{SVector{3, Float32}}();
        compobonds = Dict{Symbol, Dict{Symbol, Bondc}}()
        while !occursin("data_", rl) && !eof(compoundlibfile)
            global curcategory
            # if in cif string using " " (alt_atom_id record or ather ) change space to %% and delete " "
            if occursin("loop_", rl) 
                global curcategory
                curcategory = :absence 
            end
            split_rl = split(rl)
            # categories and attributes reord
            if length(split_rl)>0 && split_rl[1][1] == '_'
                global curcategory
                categ, atrib = split(split_rl[1], '.')
                if curcategory == :absence global curcategory = Symbol(categ) end
                if Symbol(categ) ∉ keys(categories) categories[Symbol(categ)] = Dict{Symbol, String}() end
                categories[Symbol(categ)][Symbol(atrib)] = length(split_rl) > 2 ? split_rl[2:end] : []
            end
            # if ferst substr = name of compoun it atom or bond record
            if length(split_rl)>1 && split_rl[1] == String(compoundname)
                global curcategory
                if curcategory == :_chem_comp_atom
                    if length(split_rl) != 18 
                        @error "Composition $(compoundname) is not parsing. \
                                Edit components.cif(may by in atomic record atom labeled strings with space \
                                or some adding fields) \
                                or remove pdb-structure with $(compoundname) from data"
                        error("Composition is not parsing")
                    end                    
                    atomic[Symbol(split_rl[2])] = Atomc(split_rl)
                    if !isnothing(atomic[Symbol(split_rl[2])].model_Cartn_x)
                        push!(coordcm, [atomic[Symbol(split_rl[2])].model_Cartn_x,
                                        atomic[Symbol(split_rl[2])].model_Cartn_y,
                                        atomic[Symbol(split_rl[2])].model_Cartn_z])
                    end
                    if !isnothing(atomic[Symbol(split_rl[2])].pdbx_model_Cartn_x_ideal)
                        push!(coordci, [atomic[Symbol(split_rl[2])].pdbx_model_Cartn_x_ideal,
                                        atomic[Symbol(split_rl[2])].pdbx_model_Cartn_y_ideal,
                                        atomic[Symbol(split_rl[2])].pdbx_model_Cartn_z_ideal])
                    end
                elseif curcategory == :_chem_comp_bond
                    at1 = Symbol(split_rl[2]); at2 = Symbol(split_rl[3])
                    if at1 ∉ keys(compobonds) compobonds[at1] = Dict{Symbol, Bondc}() end
                    if at2 ∉ keys(compobonds) compobonds[at2] = Dict{Symbol, Bondc}() end
                    compobonds[at1][at2] = compobonds[at2][at1] = Bondc(split_rl)
                end
            end
            rl = readline(compoundlibfile)
        end
        # bonds analitic
        # not optimal - re-calculation of link lengths, but the algorithm is simpler
        for b_atom in keys(compobonds)
            coordc = length(coordcm) > length(coordci) ? coordcm : coordci
            dists = pairwise(euclidean, [coordc[atomic[b_atom].pdbx_ordinal],],
                                         coordc[[atomic[cb_atom].pdbx_ordinal for cb_atom in keys(compobonds[b_atom])]])
            angls = pairwise(cosine_dist, 
                            coordc[[atomic[cb_atom].pdbx_ordinal for cb_atom in keys(compobonds[b_atom])]] 
                                                    .- [coordc[atomic[b_atom].pdbx_ordinal]])
            for (icb_atom, cb_atom) in enumerate(keys(compobonds[b_atom]))
                atomic[b_atom].bonds[cb_atom] = (compobonds[b_atom][cb_atom], 
                                                dists[icb_atom], 
                                                Dict(ub_atom=>angls[icb_atom, iub_atom] 
                                                    for (iub_atom, ub_atom) in enumerate(keys(compobonds[b_atom]))))
            end
        end
        # filling the target compounds set
        SetCompositions[compoundname] = (atomc = atomic, catec = categories)
    end
end

function CmpAtomic(Compound::AtomsGroup, Hflag::Bool)
    removableatoms = Set([:HA, :OXT, :H])
    global SetCompositions
    compound_atoms = Dict(ckAtom[].label_atom_id => SVector(ckAtom[].Cartn_x, ckAtom[].Cartn_y, ckAtom[].Cartn_z) 
                                     for ckAtom in Compound.childs[Atoma] if ckAtom[].type_symbol != :H || Hflag)
    if !isempty(setdiff(keys(compound_atoms), keys(SetCompositions[Compound.compoundname].atomic))) return false end
    absentatoms = setdiff(keys(SetCompositions[Compound.compoundname].atomic), atomslabels)
    if !Hflag
        for arec in deepcopy(absentatoms)
            if arec[1] == "H" delete!(absentatoms, arec) end end end
    if !isempty(setdiff(absentatoms, removableatoms)) return false end
end


end             #module StereoChem