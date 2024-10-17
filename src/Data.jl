module Data

using DataFrames, CSV
#using Glob
using Logging

using ..ProtLook: settings # settings
export settings

include("molmod/atom.jl")

#include("StereoChem.jl")
#import .StereoChem: loadCompounds, CmpAtomic

include("PDBxCIF.jl")
import .PDBxCIF: readCIF, constructMolecula
#using .AtomicProcessing

loggf = joinpath(settings[:DirOrg][:lgDir], "dlog.log")
logga = SimpleLogger(open(loggf, "w+"))

framePDB = CSV.read(
    joinpath(settings[:DirOrg][:prDir], settings[:DatLst]),
    DataFrame;
    header = true,
)
println(settings[:DatLst])
ProtsData = Dict{
    Symbol,
    @NamedTuple{
        atomic::Vector{TAtom},
        compic::Dict{Tuple{Int32,Symbol,Int32},TAtomsGroup},
        chanic::Dict{Tuple{Int32,Symbol},TPDBsChain},
        struic::Dict{Int32,TStructModel},
    }
}()
for pdbr in eachrow(framePDB)
    if Symbol(pdbr[:PDBId]) ∉ keys(ProtsData)
        print(pdbr[:PDBId], " ")
        atomica = PDBxCIF.readCIF(
            pdbr[:PDBId],
            joinpath(settings[:DirOrg][:dsDir], lowercase(pdbr[:FileName])),
            pdbr[:cif],
            pdbr[:gz],
        )
        if !ismissing(atomica)
            print(size(atomica))
            ProtsData[Symbol(pdbr[:PDBId])] =
                NamedTuple{(:atomic, :compic, :chanic, :struic)}((
                    atomica,
                    PDBxCIF.constructMolecula(atomica)...,
                ))
            #o, a, u = PDBxCIF.constructMolecula(atomica)
            #ProtsData[Symbol(pdbr[:PDBId])] = Dict(:atomic=>atomica, :compic=>o, :chanic=>a, :struic=>u)
        end
    else
        print(pdbr[:PDBId], " present ")
    end
end
println()
# print("Compounds: ")
#for pdib in keys(ProtsData)
#    for compis in keys(ProtsData[pdib].compic)
loadCompounds(
    Set([
        ProtsData[pdib].compic[compis].compoundname for pdib in keys(ProtsData) for
        compis in keys(ProtsData[pdib].compic)
    ]),
)

#framePDB = CSV.read("utils/2", DataFrame; header=true)
#for cifa in framePDB[!, :FileName]
#    PDBxCIF.readCIF(cifa, cifa, p"/home/alexandersn/WORK/DATA/assemblies/Asm/")
#end


end    #    module Data
