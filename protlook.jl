using FilePaths ; using FilePathsBase: / ; using Glob
using CodecZlib, TranscodingStreams
using CrystalInfoFramework

protmodels_path = cwd()/"protlook"
list_gf = glob("*.gz", protmodels_path)
println("numstruct ", size(list_gf))

t_single = t_all = time()
for prot_gfile in list_gf
    kinanet=Cif(read(GzipDecompressorStream(open(prot_gfile)), String), verbose=false, version=2, source=p"/home/alexandersn/WORK/julia/protlook/100d-assembly1.cif.gz")
    print(prot_gfile," ",time()-t_single, "\n")
    global t_single = time()
end
println("alltime:", time()-t_all)


