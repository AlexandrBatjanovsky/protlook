module PDBxCIF

using ..Datas: settings

using Downloads
using Logging
import TranscodingStreams: TranscodingStream as tcstream
import CodecZlib: GzipDecompressor as uzip
import OrderedCollections: OrderedDict
using DataFrames


# include("molmod/atom.jl")
# using .AtomI

"""
# 
#      downloadCif(PDBId; afname=missing, dir=missing, zip=true)
#
# загрузка cif-файла, с проверкой на наличие уже скачаной версии
# PDBId - четырёхбуквенный PDB-индетификатор
# afname - имя для сщхраняемого файла(если нет, то конструируется имя (PDBId)".cif" )
# dir - относительный путь загрузки(если нет, то используется pathtocif)
# zip - флаг архива
# возвращает труе если хорошо
"""
function downloadCif(
    urlPDBId::AbstractString,
    fname::AbstractString,
    zip::Bool = true,
)::Bool
    try
        fstream = open(fname, "w")
        Downloads.download("https://files.rcsb.org/download/" * urlPDBId, fstream)
    catch edown
        if isa(edown, SystemError)
            @warn "Cant create file $(afname)"
        elseif isa(edown, RequestError)
            @warn "Cant download file $("https://files.rcsb.org/download/"*urlPDBId)"
        else
            @warn "unhandled exception when download $(urlPDBId)"
        end
        return false
    end

    return true
end

"""

    readCIF(PDBId, afname=missing, dir=missing, zip=true)

чтение cif-файла
PDBId - четырёхбуквенный PDB-индетификатор
afname - имя для сщхраняемого файла(если нет, то конструируется имя (PDBId)".cif" )
dir - относительный путь загрузки(если нет, то используется pathtocif)
zip - флаг архива
возвращает missing или "вектор" атомов в стуктуре
"""


function readCIF(PDBId::AbstractString, afname::AbstractString, cifflag::Bool, zip::Bool)

    if !cifflag
        @warn "pdb format not supported"
        return missing
    end

    # open cif(zip or not)
    if !(isfile(afname) || downloadCif(basename(afname), afname, zip))
        @warn "missing and won't download $(PDBId)"
        return missing
    end

    try
        if zip
            global cifstream = tcstream(uzip(), open(afname, "r"))
        else
            global cifstream = open(afname, "r")
        end
    catch eread
        @warn "Cant read $(PDBId) because $(eread)"
        return missing
    end

    # cif loops split and formal examin
    # категории текущего цикла
    cur_loop_categories = OrderedDict{Symbol,OrderedDict{Symbol,Vector{String}}}()    # may be not ordered ?
    # возвращаемый список всех атомных записей в виде кортежа именнованного категориями
    # -текущего цикла
    ##atomicdata = Vector{Atoma}
    # сэт номеров циклов содержащих атомные записи (теоретически должен содержать одно значение)
    atomicloop = Dict{Int16,DataFrame}()
    numloop = 1
    for cifline in eachline(cifstream)
        if strip(cifline) == ""
            cifline = "#"
        end
        cifsplitline = split(cifline)
        if cifline == "loop_"
            numloop += 1
            empty!(cur_loop_categories)     # can be save to globalloops
        elseif cifsplitline[1][1] == '_'
            if !haskey(cur_loop_categories, Symbol(split(cifsplitline[1], ".")[1]))
                cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])] =
                    OrderedDict{Symbol,Vector{String}}()
            end
            if !haskey(
                cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])],
                Symbol(split(cifsplitline[1], ".")[2]),
            )
                cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])][Symbol(
                    split(cifsplitline[1], ".")[2],
                )] = Vector{String}(cifsplitline[2:end])
                @debug "Split loop after 2" cifsplitline[2:end]
            else
                @warn "In $(fname) categori $(split(cifsplitline[1], ".")[1]) attribute\
                       $(split(cifsplitline[1], ".")[2]) repeat"
                append!(
                    cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])][Symbol(
                        split(cifsplitline[1], ".")[2],
                    )],
                    Vector{String}(cifsplitline[2:end]),
                )
            end
        elseif cifsplitline[1] == "ATOM" || cifsplitline[1] == "HETATM"
            if !haskey(cur_loop_categories, :_atom_site)
                @warn "In $(afname) atom record with $(keys(cur_loop_categories)) categories. Text: \" $(cifline)\""
            else
                if !haskey(atomicloop, numloop)
                    if length(atomicloop) > 0
                        @warn "In $(fname) atom records in different loops"
                    end
                    atomicloop[numloop] = DataFrame(
                        [[] for i = 1:length(cur_loop_categories[:_atom_site])],
                        collect(keys(cur_loop_categories[:_atom_site]));
                        copycols = false,
                    )
                end
                push!(atomicloop[numloop], cifsplitline)
            end
        end
    end
    return (atomicloop)
end

end    #    module PDBxCIF
