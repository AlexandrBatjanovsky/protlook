module PDBxCIF

include("molmod/atom.jl")

using Distances

using ..Datas: settings
using Logging
using JLD
compoundoffset = jldopen(settings[:CmpInd], "r") do file
    read(file, "CompoundsDict")
end
SetCompositions = Dict{Symbol, @NamedTuple{atomc::Dict{Symbol, Atomc}, 
                                           bondc::Dict{Tuple{Symbol, Symbol}, (Bondc, 
                                                                               Float32, 
                                                                               Dict{Tuple{Symbol, Symbol}, Float32})}, 
                                           catec::Dict{Symbol, Dict{Symbol, Vector{Strings}}}}}()

function loadCompound(compoundname::Symbol)
    compoundlibfile = open(settings[:CmpLib], "r")
    seek(compoundlibfile, compoundoffset[compoundname])
    rl = readline(compoundlibfile)
    categories = Dict{Symbol, Dict{Symbol, Vector{Strings}}}()
    curcategory = :absence
    atomic = Dict{Symbol, Atomc}()
    bondic = Dict{Tuple{Symbol, Symbol}, Bondc}()
    while !occursin("data_", rl) && !eof(compoundlibfile)
        if occursin("loop_", rl) curcategory = :absence end
        split_rl = split(rl)
        if length(split_rl)>0 && split_rl[1] =="_"
            categ, atrib = split(split_rl[1], '.')
            if curcategory == :absence curcategory = Symbol(categ)
            if Symbol(categ) ∉ keys(categories) categories[Symbol(categ)] = Dict{Symbol, String}() end
            categories[Symbol(categ)][Symbol(atrib)] = length(split_rl) > 1 ? split_rl[2:end] : []
        end
        if length(split_rl)>1 && split_rl[1] == String(compoundname)
            if curcategory == :_chem_comp_atom
                atomic[Symbol(split_rl[end - 2])] = Atomc(split_rl)
            elseif curcategory == :_chem_comp_bond
                bondic[(Symbol(split_rl[2]), Symbol(split_rl[3]))] = (Bondc(split_rl),
                                                                      distance()
                bondic[(Symbol(split_rl[3]), Symbol(split_rl[2]))] = Bondc(split_rl)
            end
        rl = readline(compoundlibfile)
    end
    SetCompositions[compoundname] = NamedTuple{(:atomc, :bondc, :catec)}(atomic, bondic, categories)
end

function CompositionConfirmation(Compound::AtomsGroup)
    if Compound.compoundname ∉ keys(SetCompositions) loadCompound(Compound.compoundname) end
end


end             #module PDBxCIF