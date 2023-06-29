module PPDBxCIF

# using CEnum
using Dowloads
import TranscodingStreams: TranscodingStream as tcstream
import CodecZlib: GzipDecompressor as uzip
using OrderedCollections

"""
Relative path to cif-files
"""
pathtocif::AbstractPath = p"data/structs/"

"""

     downloadCif(PDBId, dir=missing, zip=true)

загрузка cif-файла, с проверкой на наличие уже скачаной версии
PDBId - четырёхбуквенный PDB-индетификатор
afname - имя для сщхраняемого файла(если нет, то конструируется имя $(PDBId)".cif" )
dir - относительный путь загрузки(если нет, то используется pathtocif)
zip - флаг архива
возвращает труе если хорошо
"""
function downloadCif(PDBId::AbstractString, afname=missing, dir::AbstractPath=missing,
                     zip::Bool=true)::Bool
    fname = PDBId*".cif"
    zip ? fname *=".gz"
    (afname === missing) && afname = fname
    if !exists(pathtocif/fname)
        try
            if (dir === missing) fstream = open(pathtocif/afname, "w")
            else fstream = open(dir/afname, "w") end
            download("https://files.rcsb.org/download/"*fname, fstream)
        catch edown
            if isa(edown, SystemError)
                @warn "Cant create file $((dir === missing) ? pathtocif/afname : dir/afname)"
            elseif isa(edown, RequestError)
                @warn "Cant download file $("https://files.rcsb.org/download/"*fname)"
            else
                @warn "unhandled exception when download $(PDBId)"
            end
            return false
        end
    end
    return true
end

"""

    readCIF(PDBId, afname=missing, dir=missing, zip=true)

чтение cif-файла
PDBId - четырёхбуквенный PDB-индетификатор
afname - имя для сщхраняемого файла(если нет, то конструируется имя $(PDBId)".cif" )
dir - относительный путь загрузки(если нет, то используется pathtocif)
zip - флаг архива
возвращает missing или "вектор" атомов в стуктуре
"""
function readCIF(PDBId::AbstractString, afname=missing, dir::AbstractPath=missing,
                 zip::Bool=true)
    # open cif(zip or not)
    fname = PDBId*".cif"
    zip ? fname *=".gz"
    (afname === missing) && afname = fname
    if !(exists(pathtocif/afname) || downloadCif(PDBId, afname, dir, zip))
        return missing
    end
    try
        if zip
            cifstream = tcstream(uzip(), open(pathtocif/afname, "r"))
        else
            open(pathtocif/afname, "r")
        end
    catch eread
        @warn "Cant read $(PDBId) because $(eread)"
        return missing
    end

    # cif loops split and formal examin
    # категории текущего цикла
    cloopcategories = OrderedDict{String, Vector{Vector{String}}}
    # возвращаемый список всех атомных записей в виде кортежа именнованного категориями
    #-текущего цикла
    atomicdata = Vector{NamedTuple}
    # сэт номеров циклов содержащих атомные записи (теоретически должен содержать одно значение)
    numatomicloop = Set{UInt8}()
    numloop = 1
    for cifline in eachline(cifstream)
        cifsplitline = split(cifline)
        if cifline == "loop_"
            numloop += 1
            empty!(currentloopcategories)
        elseif cifsplitline[1][1] == '_'
            if !haskey(cloopcategories, split(cifsplitline[1], ".")[1])
                cloopcategories[split(cifsplitline[1], ".")[1]] = Vector{Vector{String}}() end
            push!(cloopcategories[split(cifsplitline[1], ".")[1]],
                  pushfirst!(split(cifsplitline[1], ".")[2], cifsplitline[2:end]))
        elseif cifsplitline[1] == "ATOM" or cifsplitline[1] == "HETATM")
            if !haskey(cloopcategories, "_atom_site")
                    @error "In $(fname) atom record with $(keys(cloopcategories)) categories"
                return missing end
            if length(cloopcategories) != 1
                @warn "In $(fname) atom record with complex categories" end
            push!(numatomicloop, numloop)
            try push!(atomicdata, (;(Symbol.(cloopcategories["_atom_site"][:][1]) .=> cifsplitline)))
            catch enpakeatomrecort
                if isa(enpakeatomrecort, DimensionMismatch)
                    @error "In $(fname) atomic notation does not match cycle categories $(keys(cloopcategories))"
                else @error "In $(fname) atomic notation does not convert for some reason" end
                return missing
            end
        end
    end
    if length(numatomicloop)!=1 @warn "In $(fname) atom record some atomic loops" end
    return numatomicloop
end

end    #    module PDBxCIF
