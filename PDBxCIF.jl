module PDBxCIF

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
function downloadCif(urlPDBId::AbstractString, afname::AbstractString, dir::AbstractPath, zip::Bool=true)::Bool
    try
        fstream = open(dir/afname, "w")
        Downloads.download("https://files.rcsb.org/download/"*urlPDBId, fstream)
    catch edown
        if isa(edown, SystemError)
            @warn "Cant create file $(dir/afname)"
        elseif isa(edown, RequestError)
            @warn "Cant download file $("https://files.rcsb.org/download/"*urlPDBId)"
        else
            @warn "unhandled exception when download $(PDBId)"
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


function readCIF(PDBId::AbstractString, afname=nothing, dir=nothing,
                 zip::Bool=true)

    function CorrectAtomRecord(atoms_atri::OrderedDict{String, Vector{String}}, atri_flag,
                               atom_record::Vector{SubString{String}}, rettype::Type{T}, checkatris...) where T<:Union{Integer, Float32}
        results = [(indexC = ht_keyindex(atoms_atri, atri, true), valueC = atom_record[ht_keyindex(atoms_atri, atri, true)])
                   for atri in checkatris if haskey(atoms_atri, atri) && !(atom_record[ht_keyindex(atoms_atri, atri, true)] in [".", "?"])]
        if isempty(results) return nothing end
        # @debug begin
        #     if length(Set([ele.valueC for ele in results])) != 1 "Different attribute in $(PDBId) : $(atom_record), for $(checkatris)" end
        # end
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
        # @debug begin
        #     if length(Set([ele.valueC for ele in results])) != 1 "Different attribute in $(PDBId) : $(atom_record), for $(checkatris)" end
        # end
        atri_flag[results[1].indexC] = true
        return results[1].valueC
    end
    
    # open cif(zip or not)
    fname = PDBId*".cif"
    zip && (fname *=".gz")
    isnothing(afname) && (afname = fname)
    isnothing(dir) && (dir = pathtocif)
    if !(exists(dir/afname) || downloadCif(afname, afname, dir, zip))
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
    # -текущего цикла
    ##atomicdata = Vector{Atoma}
    # сэт номеров циклов содержащих атомные записи (теоретически должен содержать одно значение)
    atomicloop = Dict{Int16, @NamedTuple{categ::OrderedDict{String, OrderedDict{String, Vector{String}}}, atoms::Vector{Atoma}}}()
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
                atomicloop[numloop] = @NamedTuple{categ::OrderedDict{String, OrderedDict{String, Vector{String}}},
                                                  atoms::Vector{Atoma}}((cur_loop_categories, Vector{Atoma}())) end
            
            cur_loop_icategories = zeros(Bool, length(cur_loop_categories["_atom_site"]))
            cur_loop_icategories[1] = true
            curAtoma = Atoma(
                CorrectAtomRecord(cur_loop_categories["_atom_site"], cur_loop_icategories, cifsplitline, Int32, "id"),
                cifsplitline[1] == "HETATM",
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
                Dict()
                # Dict{String, String}(),
                # Vector{}()
            )
            push!(atomicloop[numloop].atoms, curAtoma)
        end
    end
    for loopnum in keys(atomicloop)
        global modelsdict, chainsdict, composdict
        ckmodel = atomicloop[loopnum].atoms[1].model
        modelsdict = OrderedDict{Int32, StructModel}(ckmodel=>StructModel(ckmodel,
                                                                          Dict("atoms"=>Vector{Ref{Atoma}}(),
                                                                               "compounds"=>Vector{Ref{AtomsGroup}}(),
                                                                               "chains"=>Vector{Ref{PDBsChain}}())))
        
        ckchain = (ckmodel, atomicloop[loopnum].atoms[1].chain)
        chainsdict = OrderedDict{Tuple{Int32, AbstractString}, PDBsChain}(ckchain=>
            PDBsChain(ckchain, Dict("atoms"=>Vector{Ref{Atoma}}(),
                                    "compounds"=>Vector{Ref{AtomsGroup}}()),
                      Dict("models"=>Ref{StructModel}(modelsdict[ckmodel]))))
        push!(modelsdict[ckmodel].childs["chains"], Ref(chainsdict[ckchain]))
        
        ckidcompound = (ckmodel, atomicloop[loopnum].atoms[1].chain, atomicloop[loopnum].atoms[1].idcompound)
        composdict = OrderedDict{Tuple{Int32, AbstractString, Int32}, AtomsGroup}(ckidcompound=>
            AtomsGroup(ckidcompound, atomicloop[loopnum].atoms[1].compound,
                       Dict("atoms"=>Vector{Ref{Atoma}}()),
                       Dict("chains"=>Ref{PDBsChain}(chainsdict[ckchain]),
                            "models"=>Ref{StructModel}(modelsdict[ckmodel]))))
        push!(modelsdict[ckmodel].childs["compounds"], Ref(composdict[ckidcompound]))
        push!(chainsdict[ckchain].childs["compounds"], Ref(composdict[ckidcompound]))
        
        for (ickAtomi, ckAtomi) in enumerate(atomicloop[loopnum].atoms)
            global modelsdict, chainsdict, composdict
            if ckmodel != ckAtomi.model
                if haskey(modelsdict, ckAtomi.model)
                    @warn "Model introduction $(ckmodel) in $(ckAtomi.model)"
                else modelsdict[ckAtomi.model] = StructModel(ckAtomi.model, Dict("atoms"=>Vector{Ref{Atoma}}(),
                                                                                 "compounds"=>Vector{Ref{AtomsGroup}}(),
                                                                                 "chains"=>Vector{Ref{PDBsChain}}())) end
            end
            if ckchain != (ckAtomi.model, ckAtomi.chain)
                if haskey(chainsdict, (ckAtomi.model, ckAtomi.chain))
                    @warn "Chain introduction $(ckchain) in $((ckAtomi.model, ckAtomi.chain))"
                else chainsdict[(ckAtomi.model, ckAtomi.chain)]=PDBsChain((ckAtomi.model, ckAtomi.chain),
                                                                         Dict("atoms"=>Vector{Ref{Atoma}}(),
                                                                              "compounds"=>Vector{Ref{AtomsGroup}}()),
                                                                         Dict("models"=>Ref{StructModel}(modelsdict[ckAtomi.model])))
                    push!(modelsdict[ckAtomi.model].childs["chains"], Ref(chainsdict[(ckAtomi.model, ckAtomi.chain)]))
                end
            end
            if ckidcompound != (ckAtomi.model, ckAtomi.chain, ckAtomi.idcompound)
                if haskey(composdict, (ckAtomi.model, ckAtomi.chain, ckAtomi.idcompound))
                    @warn "Compound introduction $(ckidcompound) in $(ckAtomi.idcompound)"
                else composdict[(ckAtomi.model, ckAtomi.chain, ckAtomi.idcompound)] =
                    AtomsGroup((ckAtomi.model, ckAtomi.chain,  ckAtomi.idcompound),
                               ckAtomi.compound,
                               Dict("atoms"=>Vector{Ref{Atoma}}()),
                               Dict("chains"=>Ref(chainsdict[(ckAtomi.model, ckAtomi.chain)]),
                                    "models"=>Ref(modelsdict[ckAtomi.model])))
                    push!(modelsdict[ckAtomi.model].childs["compounds"],
                          Ref(composdict[(ckAtomi.model, ckAtomi.chain, ckAtomi.idcompound)]))
                    push!(chainsdict[(ckAtomi.model, ckAtomi.chain)].childs["compounds"],
                          Ref(composdict[(ckAtomi.model, ckAtomi.chain, ckAtomi.idcompound)]))
                end
            end

            ckmodel = ckAtomi.model; ckchain = (ckAtomi.model, ckAtomi.chain)
            ckidcompound = (ckAtomi.model, ckAtomi.chain, ckAtomi.idcompound)

            push!(composdict[ckidcompound].childs["atoms"], Ref(ckAtomi))
            push!(chainsdict[ckchain].childs["atoms"], Ref(ckAtomi))
            push!(modelsdict[ckmodel].childs["atoms"], Ref(ckAtomi))

            ckAtomi.parents["model"] = Ref(modelsdict[ckmodel])
            ckAtomi.parents["chain"] = Ref(chainsdict[ckchain])
            ckAtomi.parents["compoud"] = Ref(composdict[ckidcompound])
        end
    end
    println("Models: ", length(modelsdict))
    for ckmodel in keys(modelsdict)
        println("In model ", ckmodel, " Chains ", [cha[2] for cha in keys(chainsdict) if cha[1]==ckmodel])
        for ckchain in keys(chainsdict)
            if ckchain[1] == ckmodel
                println("in chain ", ckchain, " compound ", [(comch[1], comch[2].compoundname) for comch in composdict if (comch[1][1], comch[1][2]) == ckchain])
            end
        end
    end
    println("")
    exit()
end

end    #    module PDBxCIF
