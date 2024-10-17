module PDBxCIF

export downloadCif, readCIF, constructMolecula

using ..Data: settings
using ..Data: TAtom, TStructModel, TAtomsGroup, TPDBsChain

using Downloads
import TranscodingStreams: TranscodingStream as tcstream
import CodecZlib: GzipDecompressor as uzip

using Logging
log_io = open(joinpath(settings[:DirOrg][:lgDir], "cif_log.txt"), "w+")
warn_logger = SimpleLogger(log_io)
#global_logger(warn_logger)

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
        return missing
    end
    # open cif(zip or not)
    if !(isfile(fname) || downloadCif(basename(fname), fname, zip))
        @warn "missing and won't download $(PDBId)"
        return missing
    end

    try
        if zip
            global cifstream = tcstream(uzip(), open(fname, "r"))
        else
            global cifstream = open(fname, "r")
        end
    catch eread
        @warn "Cant read $(PDBId) because $(eread)"
        return missing
    end

    # cif loops split and formal examin
    # категории текущего цикла
    cur_loop_categories = Dict{Symbol,Dict{Symbol,Vector{String}}}()
    cur_loop_atributes = Vector{Symbol}()
    # возвращаемый список всех атомных записей в виде кортежа именнованного категориями
    # -текущего цикла
    ##atomicdata = Vector{TAtom}
    # сэт номеров циклов содержащих атомные записи (теоретически должен содержать одно значение)
    atomica = Vector{TAtom}()
    atomicloop = Vector{Int16}()
    numloop = 1
    for cifline in eachline(cifstream)      # rewrite to catch read from cifstream
        if strip(cifline) == ""
            cifline = "#"
        end
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
                        Dict{Symbol,Vector{String}}()
                end
                # println(cifsplitline)
                if !haskey(
                    cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])],
                    Symbol(split(cifsplitline[1], ".")[2]),
                )
                    cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])][Symbol(
                        split(cifsplitline[1], ".")[2],
                    )] = Vector{String}(cifsplitline[2:end])
                    push!(cur_loop_atributes, Symbol(split(cifsplitline[1], ".")[2]))
                else
                    with_logger(warn_logger) do
                        @debug "In $(fname) categori $(split(cifsplitline[1], ".")[1]) attribute \
                            $(split(cifsplitline[1], ".")[2]) repeat"
                    end
                    append!(
                        cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])][Symbol(
                            split(cifsplitline[1], ".")[2],
                        )],
                        Vector{String}(cifsplitline[2:end]),
                    )
                end
            else
                with_logger(warn_logger) do
                    @warn "Loop category format not parsed $(cifline)"
                end
            end
        elseif cifsplitline[1] == "ATOM" || cifsplitline[1] == "HETATM"
            if flag_atom && !haskey(cur_loop_categories, :_atom_site)
                with_logger(warn_logger) do
                    @debug "In $(fname) atom record with $(keys(cur_loop_categories)) categories. Text: \" $(cifline)\""
                end
            else
                global flag_atom = false
                if flag_cate && !isempty(setdiff(cur_loop_atributes, fieldnames(TAtom)))
                    with_logger(warn_logger) do
                        @warn "In $(fname) Conflict names :_atom_site and TAtom: $(cur_loop_atributes))\
                            $(fieldnames(TAtom)[1:length(fieldnames(TAtom))-1])"
                    end
                    return missing
                end
                global flag_cate = false
                if numloop ∉ atomicloop
                    if length(atomicloop) > 0
                        with_logger(warn_logger) do
                            @warn "In $(fname) atom records in different loops"
                        end
                    end
                    push!(atomicloop, numloop)
                end
                if length(cifsplitline) != length(cur_loop_atributes)
                    @warn "In $(fname) bad string split of atom record: $(cifline)"
                end
                try
                    push!(
                        atomica,
                        TAtom(NamedTuple{Tuple(cur_loop_atributes)}(cifsplitline)),
                    )
                catch parserr
                    @warn "In $(fname) bad parse of atom record: $(cifline)"
                end
            end
        end
    end
    return (atomica)
end

function constructMolecula(atomiccontent::Vector{TAtom})
    global modelsdict, chainsdict, composdict
    ckmodel = atomiccontent[1].pdbx_PDB_model_num
    modelsdict = Dict{Int32,TStructModel}(
        ckmodel => TStructModel(
            ckmodel,
            Dict(
                TAtom => Vector{Ref{TAtom}}(),
                TAtomsGroup => Vector{Ref{TAtomsGroup}}(),
                TPDBsChain => Vector{Ref{TPDBsChain}}(),
            ),
        ),
    )

    ckchain = (ckmodel, atomiccontent[1].label_asym_id)
    chainsdict = Dict{Tuple{Int32,Symbol},TPDBsChain}(
        ckchain => TPDBsChain(
            ckchain,
            Dict(TAtom => Vector{Ref{TAtom}}(), TAtomsGroup => Vector{Ref{TAtomsGroup}}()),
            Dict(TStructModel => Ref{TStructModel}(modelsdict[ckmodel])),
        ),
    )
    push!(modelsdict[ckmodel].childs[TPDBsChain], Ref(chainsdict[ckchain]))

    ckcompd = (ckmodel, atomiccontent[1].label_asym_id, atomiccontent[1].label_seq_id)
    composdict = Dict{Tuple{Int32,Symbol,Int32},TAtomsGroup}(
        ckcompd => TAtomsGroup(
            ckcompd,
            atomiccontent[1].auth_comp_id,
            Dict(TAtom => Vector{Ref{TAtom}}()),
            Dict(
                TPDBsChain => Ref{TPDBsChain}(chainsdict[ckchain]),
                TStructModel => Ref{TStructModel}(modelsdict[ckmodel]),
            ),
        ),
    )
    push!(modelsdict[ckmodel].childs[TAtomsGroup], Ref(composdict[ckcompd]))
    push!(chainsdict[ckchain].childs[TAtomsGroup], Ref(composdict[ckcompd]))

    for ckAtomi in atomiccontent
        global modelsdict, chainsdict, composdict
        if ckmodel != ckAtomi.pdbx_PDB_model_num
            if haskey(modelsdict, ckAtomi.pdbx_PDB_model_num)
                with_logger(warn_logger) do
                    @warn "Model introduction $(ckAtomi.pdbx_PDB_model_num) in $(ckmodel)"
                end
            else
                global modelsdict
                modelsdict[ckAtomi.pdbx_PDB_model_num] = TStructModel(
                    ckAtomi.pdbx_PDB_model_num,
                    Dict(
                        TAtom => Vector{Ref{TAtom}}(),
                        TAtomsGroup => Vector{Ref{TAtomsGroup}}(),
                        TPDBsChain => Vector{Ref{TPDBsChain}}(),
                    ),
                )
            end
        end

        if ckchain != (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)
            if haskey(chainsdict, (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id))
                if ckAtomi.group_PDB != :HETATM
                    with_logger(warn_logger) do
                        @warn "Chain introduction $((ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)) in $(ckchain)"
                    end
                end
            else
                global chainsdict
                chainsdict[(ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)] = TPDBsChain(
                    (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id),
                    Dict(
                        TAtom => Vector{Ref{TAtom}}(),
                        TAtomsGroup => Vector{Ref{TAtomsGroup}}(),
                    ),
                    Dict(
                        TStructModel =>
                            Ref{TStructModel}(modelsdict[ckAtomi.pdbx_PDB_model_num]),
                    ),
                )
                push!(
                    modelsdict[ckAtomi.pdbx_PDB_model_num].childs[TPDBsChain],
                    Ref(chainsdict[(ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)]),
                )
            end
        end

        if ckcompd !=
           (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id, ckAtomi.label_seq_id)
            if haskey(
                composdict,
                (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id, ckAtomi.label_seq_id),
            )
                with_logger(warn_logger) do
                    @warn "Compound introduction $((ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id, ckAtomi.label_seq_id)) \
                        in $(ckcompd) "
                end
            else
                global composdict
                composdict[(
                    ckAtomi.pdbx_PDB_model_num,
                    ckAtomi.label_asym_id,
                    ckAtomi.label_seq_id,
                )] = TAtomsGroup(
                    (
                        ckAtomi.pdbx_PDB_model_num,
                        ckAtomi.label_asym_id,
                        ckAtomi.label_seq_id,
                    ),
                    ckAtomi.auth_comp_id,
                    Dict(TAtom => Vector{Ref{TAtom}}()),
                    Dict(
                        TPDBsChain => Ref(
                            chainsdict[(ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)],
                        ),
                        TStructModel => Ref(modelsdict[ckAtomi.pdbx_PDB_model_num]),
                    ),
                )
                push!(
                    modelsdict[ckAtomi.pdbx_PDB_model_num].childs[TAtomsGroup],
                    Ref(
                        composdict[(
                            ckAtomi.pdbx_PDB_model_num,
                            ckAtomi.label_asym_id,
                            ckAtomi.label_seq_id,
                        )],
                    ),
                )
                push!(
                    chainsdict[(ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)].childs[TAtomsGroup],
                    Ref(
                        composdict[(
                            ckAtomi.pdbx_PDB_model_num,
                            ckAtomi.label_asym_id,
                            ckAtomi.label_seq_id,
                        )],
                    ),
                )
            end
        end

        ckmodel = ckAtomi.pdbx_PDB_model_num
        ckchain = (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id)
        ckcompd = (ckAtomi.pdbx_PDB_model_num, ckAtomi.label_asym_id, ckAtomi.label_seq_id)

        push!(composdict[ckcompd].childs[TAtom], Ref(ckAtomi))
        push!(chainsdict[ckchain].childs[TAtom], Ref(ckAtomi))
        push!(modelsdict[ckmodel].childs[TAtom], Ref(ckAtomi))

        with_logger(warn_logger) do
            if length(ckAtomi.parents) != 0
                with_logger(warn_logger) do
                    @warn "Error atom parrent recode, $(ckAtomi)"
                end
                error("atom parrent recode")
            end
        end

        ckAtomi.parents[TStructModel] = Ref(modelsdict[ckmodel])
        ckAtomi.parents[TPDBsChain] = Ref(chainsdict[ckchain])
        ckAtomi.parents[TAtomsGroup] = Ref(composdict[ckcompd])
    end
    return (composdict, chainsdict, modelsdict)
end


end    #    module PDBxCIF
