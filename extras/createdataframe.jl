# change list PDB files to PDB Frame

using Pkg
println("Temporary environment. Please don't worry about downloaded packages")
Pkg.activate(; temp=true)
Pkg.add(["DataFrames", "CSV", "DocOpt"])

using DocOpt
argtable = """Create Struct list for ProtLook.

Usage:
  createdataframe.jl [options] file <name>
  createdataframe.jl [options] dir <path>

Options:
  -h --help           Show this screen.
  --outfile=<out>     Output file [default: outfile]
  --fpID=<fc>         Ferst position for PDBId  [default: 1].
  --spID=<sc>         Stop  position for PDBId  [default: 4].

"""
args = docopt(argtable)
println(args)
using DataFrames, CSV

if args["file"]
    framePDB = CSV.read(args["<name>"], DataFrame; header=false)
    CSV.write(args["--outfile"]*".csv",
                hcat(select(framePDB, :Column1 => (x -> SubString.(x, Ref(parse(Int64, args["--fpID"])), 
                                                                      Ref(parse(Int64, args["--spID"])))) => :PDBId),
                     rename(framePDB, 1=>:FileName)))
end