"""
module PDBCIF_File  defines loading and reading
CIF and PDB format modeles.
"""
module ProtLook


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

using JSON
#settings get/init
settings = Dict{Symbol, Union{String, Dict{Symbol, String}}}()
if !isfile(joinpath(projectDir, "settings.json"))
	println("Create zero settings")
	settings = Dict(:DirOrg=> 
						Dict(:prDir=>projectDir,
							 :scDir=>joinpath(projectDir, srcDir),
					         :wdDir=>joinpath(projectDir, dataDir),
							 :dsDir=>realpath(joinpath(projectDir, structDir)),
					         :lgDir=>realpath(joinpath(projectDir, logDir))),
					:DatLst=>
					    "PDBSelector.csv")
						
	open(joinpath(projectDir, "settings.json"), "w") do f
		JSON.print(f, settings, 2)
	end
else settings = JSON.Parser.parse(open(joinpath(projectDir, "settings.json"), "r")) end

#data init
include("Data.jl")


end    #    module ProtLook
"""
t_single = t_all = time()
for prot_gfile in list_gf
    kinanet=Cif(read(GzipDecompressorStream(open(prot_gfile)), String), verbose=false, version=2, source=p"/home/alexandersn/WORK/julia/protlook/100d-assembly1.cif.gz")
    print(prot_gfile," ",time()-t_single, "\n")
    global t_single = time()
end
println("alltime:", time()-t_all)
"""


