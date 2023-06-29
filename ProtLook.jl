"""
module PDBCIF_File  defines loading and reading
CIF and PDB format modeles.
"""
module ProtLook

using FilePaths ; using FilePathsBase: / ; using Glob
using CodecZlib, TranscodingStreams
using DataFrames, CSV
using ProgressMeter, Distributed

include("data.jl")

end    #    module ProtLook
"""
t_single = t_all = time()
for prot_gfile in list_gf
    kinanet=Cif(read(GzipDecompressorStream(open(prot_gfile)), String), verbose=false, version=2, source=p"/home/alexandersn/WORK/julia/protlook/100d-assembly1.cif.gz")
    print(prot_gfile," ",time()-t_single, "\n")
    global t_single = time()
end
println("alltime:", time()-t_all)
"""


