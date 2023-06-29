module Datas

using DataFrames, CSV
using FilePaths ; using FilePathsBase: / ; using Glob
using ProgressMeter
using Logging

include("PDBxCIF.jl")

logio = open("data_log.txt", "w")
logger = SimpleLogger(logio)

include("PPDBxCIF.jl")

framePDB = CSV.read("2", DataFrame; header=true)

pdb_download = unique(framePDB[!, :PDBId])
for pdbid in pdb_download
    println(pdbid)
    downloadCif(pdbid)
end

end    #    module Datas
