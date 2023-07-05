module PPDBxCIF

using Downloads
import TranscodingStreams: TranscodingStream as tcstream
import CodecZlib: GzipDecompressor as uzip
import OrderedCollections;
using FilePaths ; using FilePathsBase: /
    
include("molmod/atom.jl")
using .AtomI


"""
Relative path to cif-files
"""
pathtocif::AbstractPath = p"wdata/structs/"

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
function downloadCif(PDBId::AbstractString, afname=missing, dir=missing,
                     zip::Bool=true)::Bool
    fname = PDBId*".cif"
    zip && (fname *=".gz")
    (afname === missing) && (afname = fname)
    if !exists(pathtocif/fname)
        try
            if (dir === missing) fstream = open(pathtocif/afname, "w")
            else fstream = open(dir/afname, "w") end
            Downloads.download("https://files.rcsb.org/download/"*fname, fstream)
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
afname - имя для сщхраняемого файла(если нет, то конструируется имя (PDBId)".cif" )
dir - относительный путь загрузки(если нет, то используется pathtocif)
zip - флаг архива
возвращает missing или "вектор" атомов в стуктуре
"""


function readCIF(PDBId::AbstractString, afname=missing, dir=missing,
                 zip::Bool=true)

    function CorrectAtomRecord(atoms_atri::OrderedDict{String, Vector{String}}, atom_record, attributs...)::String
        # ht_keyindex not exporting, so include
        results = [atom_record[ht_keyindex(atoms_atri, atri, true)] for atri in checkatrs
                       if haskey(atoms_atri, atri) && atom_record[ht_keyindex(atoms_atri, atri, true)] not in [".", "?"]]
        if isempty(checkatrs) return missing end
    
        if length(Set(results)) != 0
            @warn "Different attribute in $(PDBId) : $(atom_record), for $(attributs)"
        end
        return results[1]
    end

    # open cif(zip or not)
    fname = PDBId*".cif"
    zip && (fname *=".gz")
    (afname === missing) && (afname = fname)
    if !(exists(pathtocif/afname) || downloadCif(PDBId, afname, dir, zip))
        return missing
    end
    try
        if zip
            global cifstream = tcstream(uzip(), open(pathtocif/afname, "r"))
        else global cifstream = open(pathtocif/afname, "r") end
    catch eread
        @warn "Cant read $(PDBId) because $(eread)"
        return missing
    end

    # cif loops split and formal examin
    # категории текущего цикла
    cloopcategories = Dict{String, OrderedDict{String, Vector{String}}}()
    # возвращаемый список всех атомных записей в виде кортежа именнованного категориями
    #-текущего цикла
    atomicdata = Vector{Atoma}
    # сэт номеров циклов содержащих атомные записи (теоретически должен содержать одно значение)
    numatomicloop = Set{UInt8}()
    numloop = 1
    for cifline in eachline(cifstream)
        cifsplitline = split(cifline)
        if cifline == "loop_"
            numloop += 1
            empty!(cloopcategories)
        elseif cifsplitline[1][1] == '_'
            if !haskey(cloopcategories, split(cifsplitline[1], ".")[1])
                cloopcategories[split(cifsplitline[1], ".")[1]] = OrderedDict{String, Vector{String}}() end
            if !haskey(cloopcategories[split(cifsplitline[1], ".")[1]], split(cifsplitline[1], ".")[2])
                cloopcategories[split(cifsplitline[1], ".")[1]][split(cifsplitline[1], ".")[2]] =
                    Vector{String}(cifsplitline[2:end])
            else
                @warn "In $(fname) categori $(split(cifsplitline[1], ".")[1]) attribute $(split(cifsplitline[1], ".")[2]) repeat"
                append!(cloopcategories[split(cifsplitline[1], ".")[1]][split(cifsplitline[1], ".")[2]], Vector{string}(cifsplitline[2:end]))
            end
        elseif cifsplitline[1] == "ATOM" || cifsplitline[1] == "HETATM"
            if !haskey(cloopcategories, "_atom_site")
                    @error "In $(fname) atom record with $(keys(cloopcategories)) categories"
                return missing end
            if length(cloopcategories) != 1
                @warn "In $(fname) atom record with complex categories" end
            push!(numatomicloop, numloop)
            # println([(str, cloopcategories["_atom_site"][collect(keys(cloopcategories["_atom_site"]))[num]]) for (num, str) in enumerate(cifsplitline)])
            println([(str, collect(keys(cloopcategories["_atom_site"]))[num]) for (num, str) in enumerate(cifsplitline)])
            
            curAtoma = Atoma(
                parse(Int32, CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "id")),
                CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "auth_atom_id", "label_atom_id"),
                CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "label_alt_id"),
                missing,
                [
                    parse(Float32, CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "Cartn_x")),
                    parse(Float32, CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "Cartn_y")),
                    parse(Float32, CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "Cartn_z"))
                ],
                parse(Float32, CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "occupancy")),
                parse(Float32, CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "B_iso_or_equiv")),
                seg!!,
                CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "type_symbol")
                parse(Float32, CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "")),
                CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "pdbx_PDB_model_num")
                CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "auth_asym_id", "label_asym_id")
                CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "auth_comp_id", "label_comp_id")
                CorrectAtomRecord(cloopcategories["_atom_site"], cifsplitline, "auth_seq_id", "label_seq_id")
                
                Dict{String, String}(),
                Vector{}()
                
            )
            # try push!(atomicdata, (;(Symbol.(cloopcategories["_atom_site"][:][1]) .=> cifsplitline)))
            # catch enpakeatomrecort
            #     if isa(enpakeatomrecort, DimensionMismatch)
            #         @error "In $(fname) atomic notation does not match cycle categories $(keys(cloopcategories))"
            #     else @error "In $(fname) atomic notation does not convert for some reason" end
            #     return missing
            # end
        end
    end
    if length(numatomicloop)!=1 @warn "In $(fname) atom record some atomic loops" end

    
    
end

end    #    module PDBxCIF
