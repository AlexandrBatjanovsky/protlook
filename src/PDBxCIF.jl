module PDBxCIF

using ..Datas: settings

using Downloads
using Logging
import TranscodingStreams: TranscodingStream as tcstream
import CodecZlib: GzipDecompressor as uzip
import OrderedCollections:OrderedDict;
#import OrderedCollections:ht_keyindex;
    
include("molmod/atom.jl")
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

    function CorrectAtomRecord(rettype::Type{T}, checkatris...) where T<:Union{Integer, Float32}
        global cur_loop_categories
        global cur_loop_icategories
        results = [findferst(==(atri), keys(cur_loop_categories[:_atom_site]) for atri in checkatris]

        returnres = nothing
        i_checkatris = 1
        while i_checkatris < length(results) + 1
            atri_value = cur_loop_categories[:_atom_site][keys(cur_loop_categories[:_atom_site])[results[i_checkatris]]]
            if !(atri_value in [".", "?"])
                returnres = tryparse(rettype, atri_value)
            if !isnothing(returnres)
                cur_loop_icategories[results[i_checkatris]] = true
                break end
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

    if !cifflag
        @warn "pdb format not supported"
        return missing end

    # open cif(zip or not)
    if !(isfile(afname) || downloadCif(basename(afname), afname, zip))
        @warn "missing and won't download $(PDBId)"
        return missing end

    try
        if zip
             global cifstream = tcstream(uzip(), open(afname, "r"))
        else global cifstream = open(afname, "r") end
    catch eread
        @warn "Cant read $(PDBId) because $(eread)"
        return missing
    end

    # cif loops split and formal examin
    # категории текущего цикла
    cur_loop_categories = OrderedDict{Symbol, OrderedDict{Symbol, Vector{String}}}()    # may be not ordered ?
    # возвращаемый список всех атомных записей в виде кортежа именнованного категориями
    # -текущего цикла
    ##atomicdata = Vector{Atoma}
    # сэт номеров циклов содержащих атомные записи (теоретически должен содержать одно значение)
    ### atomicloop = Dict{Int16, @NamedTuple{categ::OrderedDict{Symbol, OrderedDict{Symbol, Vector{String}}}, 
    # ##                                    atoms::Vector{Atoma}}}()
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
                    OrderedDict{Symbol, Vector{String}}() end
            if !haskey(cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])], 
                       Symbol(split(cifsplitline[1], ".")[2]))
                cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])]
                                   [Symdol(split(cifsplitline[1], ".")[2])] = 
                    Vector{String}(cifsplitline[2:end])
                @debug "Split loop after 2" cifsplitline[2:end]
            else
                @warn "In $(fname) categori $(split(cifsplitline[1], ".")[1]) 
                       attribute $(split(cifsplitline[1], ".")[2]) repeat"
                append!(cur_loop_categories[Symbol(split(cifsplitline[1], ".")[1])]
                                           [Symbol(split(cifsplitline[1], ".")[2])],
                        Vector{String}(cifsplitline[2:end]))
            end
        elseif cifsplitline[1] == "ATOM" || cifsplitline[1] == "HETATM"
            if !haskey(cur_loop_categories, :_atom_site)
                @error "In $(fname) atom record with $(keys(cur_loop_categories)) categories";
                return missing end
            if !haskey(atomicloop, numloop)
                if length(atomicloop) > 0 @warn "In $(fname) atom records in different loops" end
                atomicloop[numloop] = @NamedTuple{categ::OrderedDict{Symbol, OrderedDict{Symbol, Vector{String}}},
                                                  atoms::Vector{Atoma}}((cur_loop_categories, Vector{Atoma}())) end
            NamedTuple(keys(cur_loop_categories[:_atom_site]), [[] for i in 1:length(cur_loop_categories[:_atom_site])])
"""            
            curAtoma = Atoma(
                CorrectAtomRecord(cifsplitline, Int32, "id"),
                cifsplitline[1] == "HETATM",
                CorrectAtomRecord(cifsplitline, Symbol, :auth_atom_id, :label_atom_id),
                CorrectAtomRecord(cifsplitline, Symbol, :label_alt_id),
                [
                    CorrectAtomRecord(cifsplitline, Float32, :Cartn_x),
                    CorrectAtomRecord(cifsplitline, Float32, :Cartn_y),
                    CorrectAtomRecord(cifsplitline, Float32, :Cartn_z)
                ],
                CorrectAtomRecord(cur_loop_categories[:_atom_site], cur_loop_icategories,
                                  cifsplitline, Float32, :occupancy),
                CorrectAtomRecord(cur_loop_categories[:_atom_site], cur_loop_icategories,
                                  cifsplitline, Float32, :B_iso_or_equiv),
                
                CorrectAtomRecord(cur_loop_categories[:_atom_site], cur_loop_icategories,
                                  cifsplitline, Symbol, :type_symbol),
                CorrectAtomRecord(cur_loop_categories[:_atom_site], cur_loop_icategories,
                                  cifsplitline, Float32, :pdbx_formal_charge),
            
                CorrectAtomRecord(cur_loop_categories[:_atom_site], cur_loop_icategories, 
                                  cifsplitline, Int32, :pdbx_PDB_model_num),
                CorrectAtomRecord(cur_loop_categories[:_atom_site], cur_loop_icategories, 
                                  cifsplitline, Symbol, :auth_asym_id, :label_asym_id),
                CorrectAtomRecord(cur_loop_categories[:_atom_site], cur_loop_icategories, 
                                  cifsplitline, Symbol, :auth_comp_id, :label_comp_id),
                CorrectAtomRecord(cur_loop_categories[:_atom_site], cur_loop_icategories,
                                  cifsplitline, Int32, :auth_seq_id, :label_seq_id),

                cur_loop_icategories,
                Dict()
                # Dict{String, String}(),
                # Vector{}()
            )
            """
            push!(atomicloop[numloop].atoms, curAtoma)
        end
    end
    modelloop = Dict{Int16, OrderedDict{Int32, StructModel}}()
    chainloop = Dict{Int16, OrderedDict{Tuple{Int32, AbstractString}, PDBsChain}}()
    compoloop = Dict{Int16, OrderedDict{Tuple{Int32, AbstractString, Int32}, AtomsGroup}}()
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
                    if !(ckAtomi.hetatom)
                        @warn "Chain introduction $(ckchain) in $((ckAtomi.model, ckAtomi.chain))" end
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
        modelloop[loopnum] = modelsdict
        chainloop[loopnum] = chainsdict
        compoloop[loopnum] = composdict
    end
    println(PDBId)
    println("Models: ", length(modelsdict))
    for loopnum in keys(modelloop)
        for ckmodel in keys(modelloop[loopnum])
            println("In model ", ckmodel, " Chains ", [cha[2] for cha in keys(chainsdict) if cha[1]==ckmodel])
            for ckchain in keys(chainsdict)
                if ckchain[1] == ckmodel
                     # println("in chain ", ckchain, " compound ", [(comch[1], comch[2].compoundname) for comch in composdict if (comch[1][1], comch[1][2]) == ckchain])
                end
            end
        end
    end
    println("")
    return()
end

end    #    module PDBxCIF
