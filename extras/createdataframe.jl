# change list PDB files to PDB Frame

using DataFrames, CSV

framePDB = CSV.read("1", DataFrame; header=false)
CSV.write("2", hcat(select(framePDB, :Column1 => (x -> SubString.(x, Ref(1), Ref(4))) => :PDBId), rename(framePDB, 1=>:FileName)))