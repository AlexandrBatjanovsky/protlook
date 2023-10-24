module Datas

using DataFrames, CSV
using Glob
using Logging

using ..ProtLook: settings # settings
export settings

include("PDBxCIF.jl")
using .PDBxCIF
#using .AtomicProcessing

loggf = joinpath(settings[:DirOrg][:lgDir],"dlog.log")
logga = SimpleLogger(open(loggf, "w"))

framePDB = CSV.read(joinpath(settings[:DirOrg][:prDir], settings[:DatLst]), DataFrame; header=true)

ProtsData = Dict{Symbol, (Vector{Atoma}, Vector{AtomsGroup}, Vector{PDBsChain}, Vector{StructModel})}
for pdbr in eachrow(framePDB)
    print(pdbr[:PDBId], " ")
    atoma = PDBxCIF.readCIF(pdbr[:PDBId], joinpath(settings[:DirOrg][:dsDir], pdbr[:FileName]), pdbr[:cif], pdbr[:gz])
    if !ismissing(atoma) 
        model, chain, compo = PDBxCIF.constructMolecula(atoma)
    end

end

for 

#framePDB = CSV.read("utils/2", DataFrame; header=true)
#for cifa in framePDB[!, :FileName]
#    PDBxCIF.readCIF(cifa, cifa, p"/home/alexandersn/WORK/DATA/assemblies/Asm/")
#end


end    #    module Datas