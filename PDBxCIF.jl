module PPDBxCIF

using Downloads
import TranscodingStreams: TranscodingStream as tcstream
import CodecZlib: GzipDecompressor as uzip
import OrderedCollections:OrderedDict;
import OrderedCollections:ht_keyindex;
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

    function CorrectAtomRecord(atoms_atri::OrderedDict{String, Vector{String}}, atri_flag,
                               atom_record::Vector{SubString{String}}, rettype::Type{T}, checkatris...) where T<:Union{Integer, Float32}
        results = [(indexC = ht_keyindex(atoms_atri, atri, true), valueC = atom_record[ht_keyindex(atoms_atri, atri, true)])
                   for atri in checkatris if haskey(atoms_atri, atri) && !(atom_record[ht_keyindex(atoms_atri, atri, true)] in [".", "?"])]
        if isempty(results) return nothing end
        @debug begin
            if length(Set([ele.valueC for ele in results])) != 1 "Different attribute in $(PDBId) : $(atom_record), for $(checkatris)" end
        end
        returnres = nothing
        i_checkatris = 1
        while i_checkatris < length(results) + 1
            returnres = tryparse(rettype, results[i_checkatris].valueC)
            if !isnothing(returnres)
                atri_flag[results[i_checkatris].indexC] = true
                break
            end
            i_checkatris += 1
        end
        return returnres
    end

    function CorrectAtomRecord(atoms_atri::OrderedDict{String, Vector{String}}, atri_flag,
                               atom_record::Vector{SubString{String}}, rettype::Type{T}, checkatris...) where T<:AbstractString
        results = [(indexC = ht_keyindex(atoms_atri, atri, true), valueC = atom_record[ht_keyindex(atoms_atri, atri, true)])
                   for atri in checkatris if haskey(atoms_atri, atri) && !(atom_record[ht_keyindex(atoms_atri, atri, true)] in [".", "?"])]
        if isempty(results) return nothing end
        @debug begin
            if length(Set([ele.valueC for ele in results])) != 1 "Different attribute in $(PDBId) : $(atom_record), for $(checkatris)" end
        end
        atri_flag[results[1].indexC] = true
        return results[1].valueC
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
    cur_loop_categories = OrderedDict{String, OrderedDict{String, Vector{String}}}()    # may be not ordered ?
    # возвращаемый список всех атомных записей в виде кортежа именнованного категориями
    #-текущего цикла
    atomicdata = Vector{Atoma}
    # сэт номеров циклов содержащих атомные записи (теоретически должен содержать одно значение)
    atomicloop = Dict{Int16, Tuple{OrderedDict{String, OrderedDict{String, Vector{String}}}, Vector{Atoma}}}()
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
            if !haskey(cur_loop_categories, split(cifsplitline[1], ".")[1])
                cur_loop_categories[split(cifsplitline[1], ".")[1]] = OrderedDict{String, Vector{String}}() end
            if !haskey(cur_loop_categories[split(cifsplitline[1], ".")[1]], split(cifsplitline[1], ".")[2])
                cur_loop_categories[split(cifsplitline[1], ".")[1]][split(cifsplitline[1], ".")[2]] =
                    Vector{String}(cifsplitline[2:end])
            else
                @warn "In $(fname) categori $(split(cifsplitline[1], ".")[1]) attribute $(split(cifsplitline[1], ".")[2]) repeat"
                append!(cur_loop_categories[split(cifsplitline[1], ".")[1]][split(cifsplitline[1], ".")[2]],
                        Vector{String}(cifsplitline[2:end]))
            end
        elseif cifsplitline[1] == "ATOM" || cifsplitline[1] == "HETATM"
            if !haskey(cur_loop_categories, "_atom_site")
                @error "In $(fname) atom record with $(keys(cur_loop_categories)) categories";
                return missing end
            if !haskey(atomicloop, numloop)
                if length(atomicloop) > 0 @warn "In $(fname) atom records in different loops" end
                atomicloop[numloop] = (cur_loop_categories, Vector{Atoma}()) end
            
            cur_loop_icategories = zeros(Bool, length(cur_loop_categories["_atom_site"]))
            cur_loop_icategories[1] = true
            curAtoma = Atoma(
                cifsplitline[1] == "HETATM",
                CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, Int32, "id"),
                CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, String, "auth_atom_id", "label_atom_id"),
                CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, String, "label_alt_id"),
                [
                    CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, Float32, "Cartn_x"),
                    CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, Float32, "Cartn_y"),
                    CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, Float32, "Cartn_z")
                ],
                CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, Float32, "occupancy"),
                CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, Float32, "B_iso_or_equiv"),
                
                CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, String, "type_symbol"),
                CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, Float32, "pdbx_formal_charge"),
            
                CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, Int32, "pdbx_PDB_model_num"),
                CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, String, "auth_asym_id", "label_asym_id"),
                CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, String, "auth_comp_id", "label_comp_id"),
                CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, Int32, "auth_seq_id", "label_seq_id"),

                cur_loop_icategories,
                Dict("loop" => Ref(atomicloop[numloop]))
                # Dict{String, String}(),
                # Vector{}()
            )
            push!(atomicloop[numloop][2], curAtoma)
        end
    end
    for loopnum in keys(atomicloop)
        chmodel = atomicloop[loopnum][2][1].model
        modelsDict = Dict(chmodel=>sModel(chmodel, Vector{Int32}(), Dict{String, Vector{Int32}}()))
        println("Model $(chmodel)")
        for (ichAtomi, chAtomi)  in enumerate(atomicloop[loopnum][2])
            if chAtomi.model == chmodel
                push!(modelsDict[chmodel].iatoms, ichAtomi)
                chAtomi.parents["model"] = Ref(modelsDict[chmodel])
            elseif !haskey(modelsDict, chAtomi.model)
                chmodel = chAtomi.model
                println("Model $(chmodel)")
                modelsDict[chmodel] = sModel(chmodel, Vector{Int32}([ichAtomi]), Dict{String, Vector{Int32}}())
                chAtomi.parents["model"] = Ref(modelsDict[chmodel])
            else
                @warn "Model introduction $(chmodel) in $(chAtomi.model)"
                chmodel = chAtomi.model
                push!(modelsDict[chmodel].iatoms, ichAtomi)
                chAtomi.parents["model"] = Ref(modelsDict[chmodel])
            end
            if !isnothing(chAtomi.alternativ)
                if !haskey(modelsDict[chmodel].alternativ, chAtomi.alternativ)
                    modelsDict[chmodel].alternativ[chAtomi.alternativ] = Vector{Int32}()
                    println("Alternative $(chAtomi.alternativ) in model $(chmodel)")
                end
                push!(modelsDict[chmodel].alternativ[chAtomi.alternativ], ichAtomi)
            end
        end
    end        
end

end    #    module PDBxCIF
