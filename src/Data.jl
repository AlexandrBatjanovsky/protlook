module Datas

using DataFrames, CSV
using Glob
using Logging

include("PDBxCIF.jl")
using .PDBxCIF
using .AtomicProcessing

using ..ProtLook: settings
loggf = "dlog.log"
logga = SimpleLogger(open(joinpath(settings[:DirOrg][:lgDir], loggf), "w"))

framePDB = CSV.read(joinpath(settings[:DirOrg][:prDir], settings[:DatLst]), DataFrame; header=true)

pdb_download = unique(framePDB[!, :PDBId])
for pdbid in pdb_download
    println("^^^^", pdbid)
    PDBxCIF.readCIF(pdbid)
end

#framePDB = CSV.read("utils/2", DataFrame; header=true)
#for cifa in framePDB[!, :FileName]
#    PDBxCIF.readCIF(cifa, cifa, p"/home/alexandersn/WORK/DATA/assemblies/Asm/")
#end


end    #    module Datas
