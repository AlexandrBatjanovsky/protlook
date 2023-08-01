module Datas

using DataFrames, CSV
using FilePaths ; using FilePathsBase: / ; using Glob
# using ProgressMeter
using Logging

include("PDBxCIF.jl")
using .PPDBxCIF

logio = open("data_log.txt", "w")
logger = SimpleLogger(logio)

framePDB = CSV.read("utils/2pdbs_Hod.txt", DataFrame; header=true)

pdb_download = unique(framePDB[!, :PDBId])
for pdbid in pdb_download
    println("^^^^", pdbid)
    PPDBxCIF.readCIF(pdbid)
end

end    #    module Datas
