"""
module PDBCIF_File  defines loading and reading
CIF and PDB format modeles.
"""
# module PDBCIF_File

using FilePaths ; using FilePathsBase: / ; using Glob
using CodecZlib, TranscodingStreams
using DataFrames, CSV
using ProgressMeter, Distributed

@everywhere function downloadCif(PDBId::AbstractString, fname=missing, zip=true)

    (fname === missing) && (fname = PDBId)
    fname = PDBId*"."*"cif"
    if zip
        fname *=".gz"
    end
    if !ispath(pwd()*"/structs/"*fname)
        download("https://files.rcsb.org/download/"*fname, pwd()*"/structs/"*fname)
    end
end

framePDB = CSV.read("2", DataFrame; header=true)

placedir = cwd()/"structs"

if !exists(placedir)
    mkdir(placedir, recursive=true)
end
cd(placedir)

#println(typeof(unique(framePDB[!, :PDBId])[1]))
#downloadCif(unique(framePDB[!, :PDBId])[1])

pdb_download = unique(framePDB[!, :PDBId])
@sync @distributed for pdbid in pdb_download
    println(pdbid)
    downloadCif(pdbid)
end

# for pdbid_loadcif in Set(framePDB[!, :PDBId])
#     downloadCif(pdbid_loadcif)
# end



# end    # module PDBCIF_File
"""
t_single = t_all = time()
for prot_gfile in list_gf
    kinanet=Cif(read(GzipDecompressorStream(open(prot_gfile)), String), verbose=false, version=2, source=p"/home/alexandersn/WORK/julia/protlook/100d-assembly1.cif.gz")
    print(prot_gfile," ",time()-t_single, "\n")
    global t_single = time()
end
println("alltime:", time()-t_all)
"""


