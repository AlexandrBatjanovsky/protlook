"""
module PDBCIF_File  defines loading and reading
CIF and PDB format modeles.
"""
module ProtLook


export settings
#activate project envivrominent and construction
using Pkg
projectDir, srcDir = splitdir(dirname(Base.source_path()))
dataDir = "wdata"
structDir = "wdata/structs"
logDir = "wdata/logs"
if !ispath(joinpath(projectDir, structDir)) 
        mkpath(joinpath(projectDir, structDir)) end
if !ispath(joinpath(projectDir, logDir)) 
        mkpath(joinpath(projectDir, logDir)) end
Pkg.activate(projectDir)
Pkg.instantiate()

using JSON3
#settings get/init
if !isfile(joinpath(projectDir, "settings.json"))
	println("Create zero settings")
	settings = Dict(:DatLst=>
					    normpath(joinpath(projectDir, dataDir, "PDBSelector.csv")),
					:CmpLib=>
						normpath(joinpath(projectDir, dataDir, "components.cif")))

	open(joinpath(projectDir, "settings.json"), "w") do f
		JSON3.pretty(f, settings) end
else
	open(joinpath(projectDir, "settings.json"), "r") do f
		global settings = copy(JSON3.read(f)) end
end
settings[:DirOrg] = Dict(:prDir=>projectDir,
						 :scDir=>joinpath(projectDir, srcDir),
						 :wdDir=>joinpath(projectDir, dataDir),
						 :dsDir=>normpath(joinpath(projectDir, structDir)),
						 :lgDir=>normpath(joinpath(projectDir, logDir)))

#check and exit if absent list of pdb-structures for working
if !isfile(settings[:DatLst])
	@error "Need data list $(settings[:CmpLib])"
	error("Absent pdblist file")
end

using JLD
#check and download Chemical Component Dictionary and create indexing for stream
if !isfile(settings[:CmpLib])
	@info "Download absent Chemical Component Dictionary to $(settings[:CmpLib])."
	print("Download absent Chemical Component Dictionary to $(settings[:CmpLib]). Wait...")
	download("https://files.wwpdb.org/pub/pdb/data/monomers/components.cif", settings[:CmpLib])
	println("done.")
end
if !isfile(settings[:CmpLib]*".ind")
	@info "Create compounds indexing file $(settings[:CmpLib]*".ind")"
	print("Create compounds indexing file. Wait...")
	dictofcompounds = Dict{Symbol, Int64}()
	PDBlibfile = open(settings[:CmpLib], "r") 
	#datas = read(PDBlibfile, String)
	#uf = 0
	for rline in eachline(PDBlibfile)
		#println(rline, "!!", uf)
		if length(rline)>=8 && rline[1:5] == "data_"
			dictofcompounds[Symbol(rline[6:8])] = position(PDBlibfile) 
		end
	end
	jldopen(settings[:CmpLib]*".ind", "w") do indPDBlibfile
		write(indPDBlibfile, "CompoundsDict", dictofcompounds) end
	println("done.")
end

#data init
include("Data.jl")


end    #    module ProtLook