module StereoChem


export loadCompounds, CmpAtomic, CompositionConfirmation

using ..ProtLook: settings
using ..Datas: Atoma, StructModel, AtomsGroup, PDBsChain
using ..Datas: Atomc, Bondc

using Distances

using Logging
using JLD2
using StaticArrays: SVector

SetCompositions = Dict{Symbol, @NamedTuple{atomc::Dict{Symbol, Atomc},
                                           catec::Dict{Symbol, Dict{Symbol, Vector{String}}}
                                           dists::(Matrix{Float32}, Array{Float32, 3})}}()

"""
    loadCompounds(CompoundNames)

Selects from components.cif a set of elementary compounds included in the database under study (amino acid residues, heterocompounds)
CompoundNames - Set of Symbols
insert Compositions in SetCompositions(dict in module StereoChem)

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
                        @warn "Composition $(compoundname) is not parsing. \
                                Edit components.cif(may by in atomic record atom labeled strings with space \
                                or some adding fields) \
                                and pdb-structure with $(compoundname) or remove this from data"
                    end                    
                    try
                        atomic[Symbol(split_rl[2])] = Atomc(split_rl)
                    catch parseerr
                        @warn "In Compound $(compoundname) bad parse of atom record: $(rl)"
                    end
                    if !isnothing(atomic[Symbol(split_rl[2])].model_Cartn_x)
                        push!(coordcm, SVector(atomic[Symbol(split_rl[2])].model_Cartn_x,
                                               atomic[Symbol(split_rl[2])].model_Cartn_y,
                                               atomic[Symbol(split_rl[2])].model_Cartn_z))
                    end
                    if !isnothing(atomic[Symbol(split_rl[2])].pdbx_model_Cartn_x_ideal)
                        push!(coordci, SVector(atomic[Symbol(split_rl[2])].pdbx_model_Cartn_x_ideal,
                                               atomic[Symbol(split_rl[2])].pdbx_model_Cartn_y_ideal,
                                               atomic[Symbol(split_rl[2])].pdbx_model_Cartn_z_ideal))
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
            for cb_atom in keys(compobonds[b_atom])
                atomic[b_atom].bonds[cb_atom] = compobonds[b_atom][cb_atom] end
        end
        coordc = length(coordcm) > length(coordci) ? coordcm : coordci
        @debug begin 
            if length(coordc) != length(atomic) "Different num coords and atoms records in $(compoundname)" end
        end
        atomsdists = pairwise(euclidean, coordc)
        bondangles = stack([pairwise(euclidean, coordc.-[coordc[i]]) for i in length(coordc)])
        # filling the target compounds set
        SetCompositions[compoundname] = (atomc = atomic, catec = categories, dists = (atomsdists, bondangles))
    end
end

function CmpAtomic(Compound::AtomsGroup, Hflag::Bool)
    removableatoms = Set([:HA, :OXT, :H])
    absentatoms = Set{Symbol}()
    global SetCompositions
    
    ModlCompound = SetCompositions[Compound.compoundname].atomc
    ModlDistance = SetCompositions[Compound.compoundname].dists
    TestCompound = Dict{Symbol, (Int16, SVector{3, Float32})}()
    #TestCompound = Dict{Symbol, SVector{3, Float32}}()
    for aA in Compound.childs[TAtoma]
        if aA[].type_symbol != :H || Hflag
            #TestCompound[aA[].label_atom_id] = (SVector(aA[].Cartn_x, aA[].Cartn_y, aA[].Cartn_z)) end
            TestCompound[aA[].label_atom_id] = (ModlCompound[aA[].label_atom_id].pdbx_ordinal,
                                                SVector(aA[].Cartn_x, aA[].Cartn_y, aA[].Cartn_z)) end
    end
    for atomA in [aA for aA in keys(ModlCompound) if aA.type_symbol != :H || Hflag]
        if atomA in keys(TestCompound)
            atomBs = [atomB for atomB in keys(ModlCompound[atomA].bonds) if atomB in keys(TestCompound)]
            dist_diference[atomA] = [ModlDistance[ModlCompound[atomA].pdbx_ordinal, ModlCompound[atomB].pdbx_ordinal][1] 
                                                for atomB in atomBs]
                            -pairwise(euclidean, [TestCompound[atomA][2],], [TestCompound[atomB][2] for atomB in atomBs)
            angl_diference[atomA]
        else
            push!(absentatoms, atomA)
        end        
    end
end


end             #module StereoChem