"""
depricated(not need)
module progect support
- check project directories
- settings
- git
"""
module ProjectCheck_ProtLook


using Pkg
Pkg.activate(".")
emptysettingproject = Dict(
    "ProjectPathStructure" => Dict(
        "extras" => Dict(),
        "docs" => Dict(),
        "src" => Dict(),
        "wdata" => Dict("structs" => Dict()),
    ),
    "ProjectAuthor" => Dict(
        "Alexander Batjanovskii" => ("AlexandrBatjanovsky", "alexandersn@yandex.ru"),
    ),
)

isinstalled(pkg::String) =
    any(x -> x.name == pkg && x.is_direct_dep, values(Pkg.dependencies()))
if !isinstalled("JSON")
    Pkg.add("JSON")
end
using JSON
if isfile("settings.json")
    settingproject = JSON.Parser.parse(open("settings.json", "r"))
else
    settingproject = emptysettingproject
    JSON.print(open("settings.json", "w"), settingproject, 2)
end

print(keys(settingproject))

parentdir = pwd()
for i in walkdir(parentdir)
    print("PROG DEPRICATED")
end



end #ProjectCheck_ProtLook
