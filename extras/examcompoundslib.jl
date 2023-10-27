# exploring Chemical Component Dictionary file

LibCCD = Dict{AbstractString, Vector{Symbol}}()
CCDfile = open("../wdata/components.cif", "r")

atomtype = Set()

for CCDline in eachline(CCDfile)
    # println(CCDline)
    global curCategory, curComosition
    splitCCDLine = split(CCDline)
    if occursin("data_", CCDline) curComosition = CCDline[6:end] end
    if length(splitCCDLine) > 0 && occursin(".", splitCCDLine[1]) curCategory = split(splitCCDLine[1], ".")[1] end
    if length(splitCCDLine) > 1
        if splitCCDLine[1] == "_chem_comp.type" 
            global curComosition
            if join(splitCCDLine[2:end], " ") ∉ keys(LibCCD)
                LibCCD[join(splitCCDLine[2:end], " ")] = [] end
            push!(LibCCD[join(splitCCDLine[2:end], " ")], Symbol(curComosition))
        end
        if curComosition == splitCCDLine[1] && splitCCDLine[4] != "H" && curCategory == "_chem_comp_atom" 
            push!(atomtype, splitCCDLine[2])
            if splitCCDLine[2][1] == 'H' println(CCDline) end
        end
    end
end

#for cotype in keys(LibCCD)
#    println(cotype, " - ", length(LibCCD[cotype]), " : ", LibCCD[cotype]) 
#    println()
#end
for Atype in atomtype
    if Atype[1]=='H' print(":", Atype, ", ") end end