### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ ca428c92-51d8-4360-b14a-aede87919f4d
split("000 C   OXT SING N N 1")[1]

# ╔═╡ 3a2e37a6-b0aa-11ee-0e1b-83864f446fbb
begin
f=open("/home/alexandersn/Work/Proteins/protlook/wdata/components.cif", "r")
ss= ""
ucc = 0
ucl = 0
uct = 0
curcategory = ""
curcom = "%%%"
typebs = Set{Symbol}()
while ! eof(f)
	if length(ss) > 5 && ss[1:5] == "data_"
		global curcom = ss[6:end]
	end
	if ss == "loop_"
		global ucl = ucc
	end
	if length(ss) >= length("_chem_comp_bond.value_order") && ss[1:length("_chem_comp_bond.value_order")] == "_chem_comp_bond.value_order"
		global uct = ucc
	end
	if length(ss) > 1 && ss[1] == '_'
		global curcategory = split(ss, ".")[1]
	end
	if curcategory == "_chem_comp_bond" && length(ss) > 5 && split(ss)[1] == curcom
		push!(typebs, Symbol(split(ss)[uct - ucl]))
	end
	global ss=readline(f)
	global ucc = ucc + 1
end
close(f)
println(typebs)
end

# ╔═╡ Cell order:
# ╠═ca428c92-51d8-4360-b14a-aede87919f4d
# ╠═3a2e37a6-b0aa-11ee-0e1b-83864f446fbb
