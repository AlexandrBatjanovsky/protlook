module StereoChem


export loadCompounds!, CmpAtomic, CompositionConfirmation

using ..Datas: settings
using ..Datas: Atoma, StructModel, AtomsGroup, PDBsChain
using ..Datas: Atomc, Bondc

using Distances

using Logging
using JLD2
using StaticArrays: SVector
compoundoffset = load(settings[:CmpInd], "compoundoffset")

SetCompositions = Dict{Symbol, @NamedTuple{atomc::Dict{Symbol, Atomc},
                                           catec::Dict{Symbol, Dict{Symbol, Vector{String}}}}}()

function loadCompounds!(compoundnames::Set{Symbol})
    global SetCompositions
    compoundlibfile = open(settings[:CmpLib], "r")
    while !isempty(compoundnames)
        compoundname = pop!(compoundnames)
        print(compoundname, " ")
        seek(compoundlibfile, compoundoffset[compoundname])
        rl = readline(compoundlibfile)
        categories = Dict{Symbol, Dict{Symbol, Vector{String}}}()
        global curcategory = :absence
        atomic = Dict{Symbol, Atomc}()
        coordcm = Vector{SVector{3, Float32}}()
        coordci = Vector{SVector{3, Float32}}()
        compobonds = Dict{Symbol, Dict{Symbol, Bondc}}()
        while !occursin("data_", rl) && !eof(compoundlibfile)
            global curcategory
            if occursin("loop_", rl) 
                global curcategory
                curcategory = :absence 
            end
            if occursin("\"", rl)
                temprl = split(rl, "\"")
                for i in 1:length(temprl)
                    if !isodd(i)
                        temprl[i] = replace(temprl[i], " " => "%%")
                    end
                end
                rl = join(temprl)
            end
            split_rl = split(rl)
            if length(split_rl)>0 && split_rl[1][1] == '_'
                global curcategory
                categ, atrib = split(split_rl[1], '.')
                if curcategory == :absence 
                    global curcategory = Symbol(categ) 
                end
                if Symbol(categ) ∉ keys(categories) categories[Symbol(categ)] = Dict{Symbol, String}() end
                categories[Symbol(categ)][Symbol(atrib)] = length(split_rl) > 2 ? split_rl[2:end] : []
            end
            if length(split_rl)>1 && split_rl[1] == String(compoundname)
                global curcategory
                if curcategory == :_chem_comp_atom
                    println(compoundname, rl)
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

        for b_atom in keys(compobonds)
            coordc = length(coordcm) > length(coordci) ? coordcm : coordci
            dists = pairwise(euclidean, [coordc[atomic[b_atom].pdbx_ordinal],],
                                         coordc[[atomic[cb_atom].pdbx_ordinal for cb_atom in keys(compobonds[b_atom])]])
            angls = pairwise(cosine_dist, 
                            coordc[[atomic[cb_atom].pdbx_ordinal for cb_atom in keys(compobonds[b_atom])]] 
                                                    .- [coordc[atomic[b_atom].pdbx_ordinal],])
            for (icb_atom, cb_atom) in enumerate(keys(compobonds[b_atom]))
                atomic[b_atom].bonds[cb_atom] = (compobonds[b_atom][cb_atom], 
                                                dists[icb_atom], 
                                                Dict(ub_atom=>angls[icb_atom, iub_atom] 
                                                    for (iub_atom, ub_atom) in enumerate(keys(compobonds[b_atom]))))
            end
        end
        SetCompositions[compoundname] = (atomc = atomic, catec = categories)
        #print(categories)
        #println(compoundname)
    end
    # println(keys(SetCompositions))
end

function CmpAtomic(Compound::AtomsGroup, Hflag::Bool)
    removableatoms = Set([:HA, :OXT, :H])
    global SetCompositions
    atomslabels = [ckAtom[].label_atom_id for ckAtom in Compound.childs[Atoma]]
    if !isempty(setdiff(atomslabels, keys(SetCompositions[Compound.compoundname].atomic))) return false end
    absentatoms = setdiff(keys(SetCompositions[Compound.compoundname].atomic), atomslabels)
    if !Hflag
        for arec in deepcopy(absentatoms)
            if arec[1] == "H" delete!(absentatoms, arec) end end end
    if !isempty(setdiff(absentatoms, removableatoms)) return false end
end

function CompositionConfirmation(Compound::AtomsGroup)
    # "PEPTIDE LINKING" "peptide linking" "L-PEPTIDE LINKING" "L-peptide linking"
    global settings
    if Compound.compoundname ∉ keys(SetCompositions) loadCompound!(Compound.compoundname) end  # rewrite for || loadCompound!
    if CmpAtomic(Compound, settings.Hflag) && CmpStereo(Compound, settings.Hflag)
    end
end


end             #module StereoChem