ciffile = open("7spk.cif", "r") do f
    for rs in eachline(f)
        curcategory = ""
        if length(rs) > 10
            splitcategory = split(rs, ".")
            splitrs = split(rs)
            if length(splitcategory) == 2 && splitcategory[1][1] == "_" curcategory = splitcategory[1] end
            if splitrs[1] == "ATOM" || splitrs[1] == "HETATM" 
                if rs[end-1] != '?' println(rs) end
            end
        end
    end
end
