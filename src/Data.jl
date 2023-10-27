module Datas

using DataFrames, CSV
#using Glob
using Logging

using ..ProtLook: settings # settings
export settings

include("PDBxCIF.jl")
using .PDBxCIF
#using .AtomicProcessing

loggf = joinpath(settings[:DirOrg][:lgDir],"dlog.log")
logga = SimpleLogger(open(loggf, "w"))

framePDB = CSV.read(joinpath(settings[:DirOrg][:prDir], settings[:DatLst]), DataFrame; header=true)
println(settings[:DatLst])
ProtsData = Dict{Symbol, @NamedTuple{atomic::Vector{Atoma}, 
                                     compic::Dict{Tuple{Int32, Symbol, Int32}, AtomsGroup}, 
                                     chanic::Dict{Tuple{Int32, Symbol}, PDBsChain}, 
                                     struic::Dict{Int32, StructModel}}}()
for pdbr in eachrow(framePDB)
    #if Symbol(pdbr[:PDBId]) ∉ keys(ProtsData)
        print(pdbr[:PDBId], " ")
        atomica = PDBxCIF.readCIF(pdbr[:PDBId], joinpath(settings[:DirOrg][:dsDir], pdbr[:FileName]), pdbr[:cif], pdbr[:gz])
        if !ismissing(atomica) 
            print(size(atomica))
            #ProtsData[Symbol(pdbr[:PDBId])] = NamedTuple{(:atomic, :compic, :chanic, :struic)}((atomica, PDBxCIF.constructMolecula(atomica)...))
            #println(Base.summarysize(ProtsData))
        end
    #end
end

include("StereoChem.jl")

#framePDB = CSV.read("utils/2", DataFrame; header=true)
#for cifa in framePDB[!, :FileName]
#    PDBxCIF.readCIF(cifa, cifa, p"/home/alexandersn/WORK/DATA/assemblies/Asm/")
#end


end    #    module Datas