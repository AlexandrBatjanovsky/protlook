module Datas

using DataFrames, CSV
using FilePaths ; using FilePathsBase: / ; using Glob
# using ProgressMeter
using Logging

include("PDBxCIF.jl")
using .PDBxCIF

logio = open("data_log.txt", "w")
logger = SimpleLogger(logio)

framePDB = CSV.read("utils/2pdbs_Hod.txt", DataFrame; header=true)

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
