module PPDBxCIF

using CEnum
using Dowloads
using CodecZlib
"""

     downloadCif(PDBId, filename=missing, zip=true)

загрузка cif-файлов, с проверкой на наличие уже скачаной версии
PDBId - четырёхбуквенный PDB-индетификатор
zip - флаг архивной копии
"""
function downloadCif(PDBId::AbstractString, dir=missing, zip=true)
    fname = PDBId*".cif"
    if zip
        fname *=".gz"
    end
    if !ispath(pwd()*"/structs/"*fname)
        (dir === missing) ? fstream = open(fname, "w") : fstream = open(dir / fname, "w")
        download("https://files.rcsb.org/download/"*fname, fstream)
    end
end

end    #    module PDBxCIF
