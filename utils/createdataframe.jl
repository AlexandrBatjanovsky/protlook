# change list PDB files to PDB Frame
# function getpdbid(stra)
#    chop(stra, head=0, tail=17) # *"-"*chop(stra, head=13, tail=7)
# end
using DataFrames, CSV

framePDB = CSV.read("pdbs_Hod.txt", DataFrame; header=false)
CSV.write("2pdbs_Hod.txt", hcat(select(framePDB, :Column1 => (x -> SubString.(x, Ref(1), Ref(4))) => :PDBId), rename(framePDB, 1=>:FileName)))
