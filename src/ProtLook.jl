"""
module Protlook. Global. Defines and load settings and nonconstant files
runing all athers.
"""
module ProtLook

import Pkg
#Pkg.resolve()
#Pkg.instantiate()
projectDir, srcDir = splitdir(dirname(Base.source_path()))
print(projectDir, srcDir)
#Pkg.activate(projectDir)
#Pkg.instantiate()

using JSON3
using JLD2

export settings
#activate project envivrominent and construction

dataDir = "wdata"
structDir = "wdata/structs"
logDir = "wdata/logs"
if !ispath(joinpath(projectDir, structDir))
    mkpath(joinpath(projectDir, structDir))
end
if !ispath(joinpath(projectDir, logDir))
    mkpath(joinpath(projectDir, logDir))
end

#settings get/init
if !isfile(joinpath(projectDir, "settings.json"))
    println("Create zero settings")
    settings = Dict(
        :DatLst => normpath(joinpath(projectDir, dataDir, "outfile.csv")),    #PDBSelector_f3b70
        :CmpLib => normpath(joinpath(projectDir, dataDir, "components.cif")),
        :CmpInd => normpath(joinpath(projectDir, dataDir, "components.cif.jld2")),
        :DirOrg => Dict(
            :prDir => projectDir,
            :scDir => joinpath(projectDir, srcDir),
            :wdDir => joinpath(projectDir, dataDir),
            :dsDir => normpath(joinpath(projectDir, structDir)),
            :lgDir => normpath(joinpath(projectDir, logDir)),
        ),
    )

    open(joinpath(projectDir, "settings.json"), "w") do f
        JSON3.pretty(f, settings)
    end
else
    open(joinpath(projectDir, "settings.json"), "r") do f
        global settings = copy(JSON3.read(f))
    end
end
settings[:DirOrg] = Dict(
    :prDir => projectDir,
    :scDir => joinpath(projectDir, srcDir),
    :wdDir => joinpath(projectDir, dataDir),
    :dsDir => normpath(joinpath(projectDir, structDir)),
    :lgDir => normpath(joinpath(projectDir, logDir)),
)

#check and exit if absent list of pdb-structures for working
if !isfile(settings[:DatLst])
    @error "Need data list $(settings[:DatLst])"
    error("Absent pdblist file")
end

#check and download Chemical Component Dictionary and create indexing for stream
if !isfile(settings[:CmpLib])
    @info "Download absent Chemical Component Dictionary to $(settings[:CmpLib])."
    print("Download absent Chemical Component Dictionary to $(settings[:CmpLib]). Wait...")
    download(
        "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif",
        settings[:CmpLib],
    )
    println("done.")
end

# create seek(read) indexes of compouds in Chemical Component Dictionary and save dict of its to JLD2
if !isfile(settings[:CmpInd])
    @info "Create compounds indexing file $(settings[:CmpLib]*".ind")"
    print("Create compounds indexing file. Wait...")
    dictofcompounds = Dict{Symbol,Int64}()
    PDBlibfile = open(settings[:CmpLib], "r")
    for rline in eachline(PDBlibfile)
        if length(rline) >= 5 && rline[1:5] == "data_"
            dictofcompounds[Symbol(strip(rline[6:end], [' ', '\n']))] = position(PDBlibfile)
        end
    end
    jldopen(settings[:CmpInd], "w"; compress = true) do indPDBlibfile
        indPDBlibfile["compoundoffset"] = dictofcompounds
    end
    println("done.")
end

#data init
include("Data.jl")

end    #    module ProtLook
