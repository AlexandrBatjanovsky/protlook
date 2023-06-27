module Datas

using DataFrames, CSV
using FilePaths ; using FilePathsBase: / ; using Glob
using TranscodingStreams, CodecZlib
using ProgressMeter

include("PPDBxCIF.jl")

framePDB = CSV.read("2", DataFrame; header=true)

placedir = cwd()/"structs"

if !exists(placedir)
    mkdir(placedir, recursive=true)
end
cd(placedir)

pdb_download = unique(framePDB[!, :PDBId])
@sync @distributed for pdbid in pdb_download
    println(pdbid)
    downloadCif(pdbid)
end





end    #    module Datas
