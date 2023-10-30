module PDBxCIF

using ..Datas: settings

using Downloads
import TranscodingStreams: TranscodingStream as tcstream
import CodecZlib: GzipDecompressor as uzip

using Logging
log_io = open("logfile.txt", "w")
debug_logger = ConsoleLogger(log_io, Logging.Debug)
#global_logger(debug_logger)

include("molmod/atom.jl")

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
function downloadCif(urlPDBId::AbstractString, fname::AbstractString, zip::Bool=true)::Bool
    try
        fstream = open(fname, "w")
        Downloads.download("https://files.rcsb.org/download/"*urlPDBId, fstream)
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

    print(" download ")
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
function readCIF(PDBId::AbstractString, fname::AbstractString, cifflag::Bool, zip::Bool)
                 
    if !cifflag
        @warn "pdb format not supported"
        return missing end
    # open cif(zip or not)
    if !(isfile(fname) || downloadCif(basename(fname), fname, zip))
        @warn "missing and won't download $(PDBId)"
        return missing end

    try
        if zip
            global cifstream = tcstream(uzip(), open(fname, "r"))
        else global cifstream = open(fname, "r") end
    catch eread
        @warn "Cant read $(PDBId) because $(eread)"
        return missing
    end

    # cif loops split and formal examin
    # категории текущего цикла
    cur_loop_categories = Dict{Symbol, Dict{Symbol, Vector{String}}}()
    cur_loop_atributes = Vector{Symbol}()
    # возвращаемый список всех атомных записей в виде кортежа именнованного категориями
    # -текущего цикла
    ##atomicdata = Vector{Atoma}
    # сэт номеров циклов содержащих атомные записи (теоретически должен содержать одно значение)
    atomica = Vector{Atoma}()
    atomicloop = Vector{Int16}()
    numloop = 1
    for cifline in eachline(cifstream)      # rewrite to catch read from cifstream
        if strip(cifline) == "" cifline = "#" end
        cifsplitline = split(cifline)
        if cifline == "loop_"
            global flag_atom = true
            global flag_cate = true 
            numloop += 1
            empty!(cur_loop_categories)     # can be save to globalloops
            empty!(cur_loop_atributes)
        elseif cifsplitline[1][1] == '_' 
            if '.' ∈ cifsplitline[1]
                if !haskey(cur_loop_categories, Symbol(split(cifsplitline[1], ".")[1]))
                    cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])] = 
                        Dict{Symbol, Vector{String}}() end
                # println(cifsplitline)
                if !haskey(cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])],
                        Symbol(split(cifsplitline[1], ".")[2]))
                    cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])][
                                        Symbol(split(cifsplitline[1], ".")[2])] =
                        Vector{String}(cifsplitline[2:end])
                    push!(cur_loop_atributes, Symbol(split(cifsplitline[1], ".")[2]))
                else
                    @warn "In $(fname) categori $(split(cifsplitline[1], ".")[1]) attribute\
                        $(split(cifsplitline[1], ".")[2]) repeat"
                    append!(cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])][
                                                Symbol(split(cifsplitline[1], ".")[2])],
                            Vector{String}(cifsplitline[2:end]))
                end
            else
                @warn "Loop category format not parsed $(cifline)"
            end
        elseif cifsplitline[1] == "ATOM" || cifsplitline[1] == "HETATM"
            if flag_atom && !haskey(cur_loop_categories, :_atom_site)
                @warn "In $(fname) atom record with $(keys(cur_loop_categories)) categories. Text: \" $(cifline)\""
            else
                global flag_atom = false
                if flag_cate && !isempty(setdiff(cur_loop_atributes, fieldnames(Atoma)))
                    @warn "In $(fname) Conflict names :_atom_site and Atoma: $(cur_loop_atributes))\
                    $(fieldnames(Atoma)[1:length(fieldnames(Atoma))-1])"
                    return missing
                end
                global flag_cate = false
                if numloop ∉ atomicloop
                    if length(atomicloop) > 0 @warn "In $(fname) atom records in different loops" end
                    push!(atomicloop, numloop)
                end
                push!(atomica, Atoma(NamedTuple{Tuple(cur_loop_atributes)}(cifsplitline)))
            end
        end
    end
    return(atomica)

end

function constructMolecula(atomiccontent::Vector{Atoma})
    global modelsdict, chainsdict, composdict
    ckmodel = atomiccontent[1].pdbx_PDB_model_num
    modelsdict = Dict{Int32, StructModel}(ckmodel=>
                                          StructModel(ckmodel,
                                                      Dict(Atoma=>Vector{Ref{Atoma}}(),
                                                           AtomsGroup=>Vector{Ref{AtomsGroup}}(),
                                                           PDBsChain=>Vector{Ref{PDBsChain}}())))
    
    ckchain = (ckmodel, atomiccontent[1].label_asym_id)
    chainsdict = Dict{Tuple{Int32, Symbol}, PDBsChain}(ckchain=>
                                                       PDBsChain(ckchain, 
                                                                 Dict(Atoma=>Vector{Ref{Atoma}}(),
                                                                      AtomsGroup=>Vector{Ref{AtomsGroup}}()),
                                                                 Dict(StructModel=>Ref{StructModel}(modelsdict[ckmodel]))))
    push!(modelsdict[ckmodel].childs[PDBsChain], Ref(chainsdict[ckchain]))

    ckcompd = (ckmodel, atomiccontent[1].label_asym_id, atomiccontent[1].label_seq_id)
    composdict = Dict{Tuple{Int32, Symbol, Int32}, AtomsGroup}(ckcompd=>
        AtomsGroup(ckcompd, atomiccontent[1].auth_comp_id,
                   Dict(Atoma=>Vector{Ref{Atoma}}()),
                   Dict(PDBsChain=>Ref{PDBsChain}(chainsdict[ckchain]),
                        StructModel=>Ref{StructModel}(modelsdict[ckmodel]))))
    push!(modelsdict[ckmodel].childs[AtomsGroup], Ref(composdict[ckcompd]))
    push!(chainsdict[ckchain].childs[AtomsGroup], Ref(composdict[ckcompd]))

    for ckAtomi in atomiccontent
        global modelsdict, chainsdict, composdict
        if ckmodel != ckAtomi.pdbx_PDB_model_num
            if haskey(modelsdict, ckAtomi.pdbx_PDB_model_num)
                @warn "Model introduction $(ckAtomi.pdbx_PDB_model_num) in $(ckmodel)"
            else
                global modelsdict
                modelsdict[ckAtomi.pdbx_PDB_model_num] = StructModel(ckAtomi.pdbx_PDB_model_num,
                                                                     Dict(Atoma=>Vector{Ref{Atoma}}(),
                                                                          AtomsGroup=>Vector{Ref{AtomsGroup}}(),
                                                                          PDBsChain=>Vector{Ref{PDBsChain}}())) 
            end
        end

        if ckchain != (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)
            if haskey(chainsdict, (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id))
                if ckAtomi.group_PDB != :HETATM
                    @warn "Chain introduction $((ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)) in $(ckchain)" end
            else 
                global chainsdict
                chainsdict[(ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)]=PDBsChain(
                    (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id),
                    Dict(Atoma=>Vector{Ref{Atoma}}(),
                         AtomsGroup=>Vector{Ref{AtomsGroup}}()),
                    Dict(StructModel=>Ref{StructModel}(modelsdict[ckAtomi.pdbx_PDB_model_num])))
                push!(modelsdict[ckAtomi.pdbx_PDB_model_num].childs[PDBsChain], 
                      Ref(chainsdict[(ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)]))
            end
        end

        if ckcompd != (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id, ckAtomi.label_seq_id)
            if haskey(composdict, (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id, ckAtomi.label_seq_id))
                print(ckAtomi)
                @warn "Compound introduction $((ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id, ckAtomi.label_seq_id))\
                             in $(ckcompd) "
            else
                global composdict 
                composdict[(ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id, ckAtomi.label_seq_id)] =
                AtomsGroup((ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id, ckAtomi.label_seq_id),
                           ckAtomi.auth_comp_id,
                           Dict(Atoma=>Vector{Ref{Atoma}}()),
                           Dict(PDBsChain=>Ref(chainsdict[(ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)]),
                                StructModel=>Ref(modelsdict[ckAtomi.pdbx_PDB_model_num])))
                push!(modelsdict[ckAtomi.pdbx_PDB_model_num].childs[AtomsGroup],
                      Ref(composdict[(ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id, ckAtomi.label_seq_id)]))
                push!(chainsdict[(ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)].childs[AtomsGroup],
                      Ref(composdict[(ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id, ckAtomi.label_seq_id)]))
            end
        end

        ckmodel = ckAtomi.pdbx_PDB_model_num; ckchain = (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)
        ckcompd = (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id, ckAtomi.label_seq_id)
        
        push!(composdict[ckcompd].childs[Atoma], Ref(ckAtomi))
        push!(chainsdict[ckchain].childs[Atoma], Ref(ckAtomi))
        push!(modelsdict[ckmodel].childs[Atoma], Ref(ckAtomi))
        
        with_logger(debug_logger) do
            if length(ckAtomi.parents)!=0 
                @debug begin
                    println("Error atom parrent recode", ckAtomi)
                    error("atom parrent recode")
                end
            end
        end 

        ckAtomi.parents[StructModel] = Ref(modelsdict[ckmodel])
        ckAtomi.parents[PDBsChain] = Ref(chainsdict[ckchain])
        ckAtomi.parents[AtomsGroup] = Ref(composdict[ckcompd])
    end
    return (composdict, chainsdict, modelsdict)
end

    
end    #    module PDBxCIF
