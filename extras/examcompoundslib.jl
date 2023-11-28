# exploring Chemical Component Dictionary file

LibCCD = Dict{Symbol, Set{Symbol}}(:H=>Set())
CCDfile = open("../wdata/components.cif", "r")

curCategory = ""
CCDline = ""
for CCDline in eachline(CCDfile)
    #println(CCDline)
    global curCategory, curComosition
    if length(CCDline) > 1
        if startswith(CCDline, "data_") 
            curComosition = CCDline[6:end] end
        if length(CCDline) > 0 && CCDline[1] == '_'
            curCategory = split(CCDline, ".")[1] 
        end
        if startswith(CCDline, curComosition) && curCategory == "_chem_comp_atom"
            if length(split(CCDline))!=18 println(length(split(CCDline)), " ", CCDline)
            elseif split(CCDline)[4] == "H"
                push!(LibCCD[:H], Symbol(split(CCDline)[2]))
            end
        end
    end
end

println(LibCCD[:H])