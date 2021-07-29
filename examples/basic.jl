### A Pluto.jl notebook ###
# v0.15.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 57d10945-b2b3-4c2f-8ae2-16db9a488ee6
begin
	# Select dropdowns
	using PlutoUI
	
	using PowerModelsDistribution 
	
	# Solvers
	using Ipopt
	using Cbc
	using SCS
	using Juniper
end

# ╔═╡ f9e56b1c-a6c4-4704-8d37-151caa33a1cc
md"""
# Using PowerModelsDistribution
"""

# ╔═╡ 00daa957-0a74-4609-9a1b-ea2b788e45a3
md"""

## Select a case file:

$(@bind case_file Select(
	[joinpath(dirname(pathof(PowerModelsDistribution)), "../test/data/opendss/case3_balanced.dss") => "case3_balanced",
	joinpath(dirname(pathof(PowerModelsDistribution)), "../test/data/opendss/case3_unbalanced.dss") => "case3_unbalanced",
	joinpath(dirname(pathof(PowerModelsDistribution)), "../test/data/opendss/ut_trans_2w_yy.dss") => "ut_trans_2w_yy",
	]
))
"""

# ╔═╡ e8d28cbb-69fb-44de-81f5-86dc903df229
md"""
## Parse case $(split(basename(case_file), ".")[1]):
"""

# ╔═╡ ab88f189-3093-4a54-80ef-9120e9a20202
data = parse_file(case_file)

# ╔═╡ 79e5480a-aecb-49fb-ba7e-d1904a13b58d
begin
	forms_for_case = Dict(
		"case3_balanced.dss" => ["AC-rectangular", "AC-polar", "Current-Voltage", "LinDistFlow", "DC", "Network Flow Approximation (NFA)", "Second Order Cone (SOC)", "Semi-definite (SDP)"],
		"case3_unbalanced.dss" => ["AC-rectangular", "AC-polar", "Current-Voltage", "LinDistFlow", "DC", "Network Flow Approximation (NFA)", "Second Order Cone (SOC)", "Semi-definite (SDP)"],
		"ut_trans_2w_yy.dss" => ["AC-rectangular", "AC-polar", "Current-Voltage", "LinDistFlow"]
		)[basename(case_file)]
	md"""
## Select a formulation to use with $(split(basename(case_file), ".")[1]):

$(@bind form_name Select(forms_for_case))
"""
end

# ╔═╡ c139b036-c4f1-4a1e-aed1-5fc77e2ab00c
formulation = Dict(
	"AC-rectangular" => ACRUPowerModel, 
	"AC-polar"=>ACPUPowerModel, 
	"Current-Voltage"=>IVRUPowerModel, 
	"LinDistFlow"=>LPUBFDiagPowerModel, 
	"DC" => DCPUPowerModel,
	"Network Flow Approximation (NFA)" => NFAUPowerModel,
	"Second Order Cone (SOC)" => SOCConicUBFPowerModel,
	"Semi-definite (SDP)" => SDPUBFPowerModel,
)[form_name]

# ╔═╡ 67033525-d437-4d12-bfe9-c7acb6df1bc2
begin
	problems_for_form = Dict(
		ACPUPowerModel => ["Optimal Power Flow (OPF)", "Power Flow (PF)", "Maximal Load Delivery (MLD)"],
		ACRUPowerModel => ["Optimal Power Flow (OPF)", "Power Flow (PF)", "Maximal Load Delivery (MLD)"],
		IVRUPowerModel => ["Optimal Power Flow (OPF)", "Power Flow (PF)"],
		LPUBFDiagPowerModel => ["Optimal Power Flow (OPF)", "Power Flow (PF)", "Maximal Load Delivery (MLD)"],
		NFAUPowerModel => ["Optimal Power Flow (OPF)", "Power Flow (PF)", "Maximal Load Delivery (MLD)"],
		DCPUPowerModel => ["Optimal Power Flow (OPF)", "Power Flow (PF)", "Maximal Load Delivery (MLD)"],
		SOCConicUBFPowerModel => ["Optimal Power Flow (OPF)", "Power Flow (PF)"],
		SDPUBFPowerModel => ["Optimal Power Flow (OPF)", "Power Flow (PF)"],

		)[formulation]

	md"""
## Select a problem to solve using $(form_name):

$(@bind prob_name Select(problems_for_form))
"""
end

# ╔═╡ 2d0f68f9-c336-4120-85bb-e781b55bcb5f
problem = Dict(
	"Optimal Power Flow (OPF)" => solve_mc_opf, 
	"Power Flow (PF)" => solve_mc_pf, 
	"Maximal Load Delivery (MLD)" => solve_mc_mld,
)[prob_name]

# ╔═╡ d8749fe0-3915-4b63-a9e3-92726b590962
begin
	solvers_for_form_problem = Dict(
		(ACPUPowerModel,solve_mc_opf) => ["Ipopt", "Juniper"],
		(ACPUPowerModel,solve_mc_pf) => ["Ipopt", "Juniper"],
		(ACPUPowerModel,solve_mc_mld) => ["Ipopt", "Juniper"],
		
		(ACRUPowerModel,solve_mc_opf) => ["Ipopt", "Juniper"],
		(ACRUPowerModel,solve_mc_pf) => ["Ipopt", "Juniper"],
		(ACRUPowerModel,solve_mc_mld) => ["Ipopt", "Juniper"],
		
		(IVRUPowerModel,solve_mc_opf) => ["Ipopt", "Juniper"],
		(IVRUPowerModel,solve_mc_pf) => ["Ipopt", "Juniper"],
		(IVRUPowerModel,solve_mc_mld) => ["Ipopt", "Juniper"],
		
		(LPUBFDiagPowerModel,solve_mc_opf) => ["Juniper", "Ipopt"],
		(LPUBFDiagPowerModel,solve_mc_pf) => ["Juniper", "Ipopt"],
		(LPUBFDiagPowerModel,solve_mc_mld) => ["Juniper", "Ipopt"],
		
		(NFAUPowerModel,solve_mc_opf) => ["Cbc", "Ipopt", "Juniper"],
		(NFAUPowerModel,solve_mc_pf) => ["Cbc", "Ipopt", "Juniper"],
		(NFAUPowerModel,solve_mc_mld) => ["Ipopt", "Juniper"],

		(DCPUPowerModel,solve_mc_opf) => ["Ipopt", "Juniper"],
		(DCPUPowerModel,solve_mc_pf) => ["Cbc", "Ipopt", "Juniper"],
		(DCPUPowerModel,solve_mc_mld) => ["Ipopt", "Juniper"],
		
		(SOCConicUBFPowerModel,solve_mc_opf) => ["SCS"],
		(SOCConicUBFPowerModel,solve_mc_pf) => ["SCS"],
		
		(SDPUBFPowerModel,solve_mc_opf) => ["SCS"],
		(SDPUBFPowerModel,solve_mc_pf) => ["SCS"],
		
	)[(formulation,problem)]
	
	md"""
## Select a solver to solve $(prob_name) with $(form_name):

$(@bind solver_name Select(solvers_for_form_problem))
"""
end

# ╔═╡ f22095fa-0be0-4b78-86f7-628a015f91f8
solver = Dict(
	"Ipopt" => optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0),
	"Cbc" => optimizer_with_attributes(Cbc.Optimizer, "logLevel"=>0),
	"SCS" => optimizer_with_attributes(SCS.Optimizer, "max_iters"=>20000, "eps"=>1e-5, "alpha"=>0.4, "verbose"=>0),
	"Juniper" => optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-4, "print_level"=>0), "mip_solver"=>optimizer_with_attributes(Cbc.Optimizer, "logLevel"=>0), "log_levels"=>[])
)[solver_name]

# ╔═╡ 8f4fd610-d4b7-45bf-9473-65f7720cb224
md"""
## Run $(prob_name)
"""

# ╔═╡ 01502fce-b6b9-4696-8562-8f674f1575f8
result = problem(data, formulation, solver)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[compat]
Cbc = "~0.8.0"
Ipopt = "~0.7.0"
Juniper = "~0.7.0"
PlutoUI = "~0.7.9"
PowerModelsDistribution = "~0.11.4"
SCS = "~0.7.1"

[deps]
Cbc = "9961bab8-2fa3-5c5a-9d89-47fab24efd76"
Ipopt = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
Juniper = "2ddba703-00a4-53a7-87a5-e8b9971dde84"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PowerModelsDistribution = "d7431456-977f-11e9-2de3-97ff7677985e"
SCS = "c946c3f1-0d1f-5ce8-9dea-7daa1f7e2d13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ASL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "370cafc70604b2522f2c7cf9915ebcd17b4cd38b"
uuid = "ae81ac8f-d209-56e5-92de-9978fef736f9"
version = "0.1.2+0"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Statistics", "UUIDs"]
git-tree-sha1 = "01ca3823217f474243cc2c8e6e1d1f45956fe872"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.0.0"

[[BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[Cbc]]
deps = ["BinaryProvider", "CEnum", "Cbc_jll", "Libdl", "MathOptInterface", "SparseArrays"]
git-tree-sha1 = "8da071e2f96de95decbd551d54db9ac82c0bf4f6"
uuid = "9961bab8-2fa3-5c5a-9d89-47fab24efd76"
version = "0.8.0"

[[Cbc_jll]]
deps = ["Artifacts", "Cgl_jll", "Clp_jll", "CoinUtils_jll", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS32_jll", "Osi_jll", "Pkg"]
git-tree-sha1 = "5d24ea46edebb4f067b180cfc092a45d4d2c48b3"
uuid = "38041ee0-ae04-5750-a4d2-bb4d0d83d27d"
version = "2.10.5+2"

[[Cgl_jll]]
deps = ["Artifacts", "Clp_jll", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0c41f887cca7e95e7ca1cda5bd6e82e6f43b317d"
uuid = "3830e938-1dd0-5f3e-8b8e-b3ee43226782"
version = "0.60.2+6"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "be770c08881f7bb928dfd86d1ba83798f76cf62a"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "0.10.9"

[[Clp_jll]]
deps = ["Artifacts", "CoinUtils_jll", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS32_jll", "Osi_jll", "Pkg"]
git-tree-sha1 = "d9eca9fa2435959b5542b13409a8ec5f64c947c8"
uuid = "06985876-5285-5a41-9fcb-8948a742cc53"
version = "1.17.6+7"

[[CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[CoinUtils_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "5186155a8609b71eae7e104fa2b8fbf6ecd5d9bb"
uuid = "be027038-0da8-5614-b30d-e42594cb92df"
version = "2.11.3+4"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dc7dedc2c2aa9faf59a55c622760a25cbefbe941"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.31.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4437b64df1e0adccc3e5d1adbc3ac741095e4677"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.9"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "214c3fcac57755cfda163d91c58893a8723f93e9"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.0.2"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[ExprTools]]
git-tree-sha1 = "10407a39b87f29d47ebaca8edbc75d7c302ff93e"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.3"

[[EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "0fa3b52a04a4e210aeb1626def9c90df3ae65268"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.1.0"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "e2af66012e08966366a43251e1fd421522908be6"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.18"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "86ed84701fbfd1142c9786f8e53c595ff5a4def9"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.10"

[[InfrastructureModels]]
deps = ["JuMP", "MathOptInterface", "Memento"]
git-tree-sha1 = "14ad999045d7ea263cb909b1942814850227d457"
uuid = "2030c09a-7f63-5d83-885d-db604e0e9cc0"
version = "0.6.1"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Ipopt]]
deps = ["BinaryProvider", "Ipopt_jll", "Libdl", "LinearAlgebra", "MathOptInterface", "MathProgBase"]
git-tree-sha1 = "380786b4929b8d18d76e909c6b2eca355b7c3bd6"
uuid = "b6b21f68-93f8-5de0-b562-5493be1d77c9"
version = "0.7.0"

[[Ipopt_jll]]
deps = ["ASL_jll", "Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "MUMPS_seq_jll", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "82124f27743f2802c23fcb05febc517d0b15d86e"
uuid = "9cc047cb-c261-5740-88fc-0cf96f7bdcc7"
version = "3.13.4+2"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[JSONSchema]]
deps = ["HTTP", "JSON", "ZipFile"]
git-tree-sha1 = "b84ab8139afde82c7c65ba2b792fe12e01dd7307"
uuid = "7d188eb4-7ad8-530c-ae41-71a32a6d4692"
version = "0.3.3"

[[JuMP]]
deps = ["Calculus", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MathOptInterface", "MutableArithmetics", "NaNMath", "Printf", "Random", "SparseArrays", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "8dfc5df8aad9f2cfebc8371b69700efd02060827"
uuid = "4076af6c-e467-56ae-b986-b466b2749572"
version = "0.21.8"

[[Juniper]]
deps = ["Distributed", "JSON", "LinearAlgebra", "MathOptInterface", "MutableArithmetics", "Random", "Statistics"]
git-tree-sha1 = "6fa2bc3df512f21c85aa82c403fa23446a1deb4e"
uuid = "2ddba703-00a4-53a7-87a5-e8b9971dde84"
version = "0.7.0"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "1ba664552f1ef15325e68dc4c05c3ef8c2d5d885"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.2.4"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "dfeda1c1130990428720de0024d4516b1902ce98"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.7"

[[METIS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "2dc1a9fc87e57e32b1fc186db78811157b30c118"
uuid = "d00139f3-1899-568f-a2f0-47f597d42d70"
version = "5.1.0+5"

[[MUMPS_seq_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "METIS_jll", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "1a11a84b2af5feb5a62a820574804056cdc59c39"
uuid = "d7ed1dd3-d0ae-5e8e-bfb4-87a502085b8d"
version = "5.2.1+4"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6a8a2a625ab0dea913aba95c11370589e0239ff0"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "JSON", "JSONSchema", "LinearAlgebra", "MutableArithmetics", "OrderedCollections", "SparseArrays", "Test", "Unicode"]
git-tree-sha1 = "575644e3c05b258250bb599e57cf73bbf1062901"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "0.9.22"

[[MathProgBase]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9abbe463a1e9fc507f12a69e7f29346c2cdc472c"
uuid = "fdba3010-5040-5b88-9595-932c9decdf73"
version = "0.7.8"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Memento]]
deps = ["Dates", "Distributed", "JSON", "Serialization", "Sockets", "Syslogs", "Test", "TimeZones", "UUIDs"]
git-tree-sha1 = "19650888f97362a2ae6c84f0f5f6cda84c30ac38"
uuid = "f28f55f0-a522-5efc-85c2-fe41dfb9b2d9"
version = "1.2.0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "916b850daad0d46b8c71f65f719c49957e9513ed"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.1"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3927848ccebcc165952dc0d9ac9aa274a87bfe01"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.20"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OpenBLAS32_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ba4a8f683303c9082e84afba96f25af3c7fb2436"
uuid = "656ef2d0-ae68-5445-9ca0-591084a874a2"
version = "0.3.12+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Osi_jll]]
deps = ["Artifacts", "CoinUtils_jll", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS32_jll", "Pkg"]
git-tree-sha1 = "ef540e28c9b82cb879e33c0885e1bbc9a1e6c571"
uuid = "7da25872-d9ce-5375-a4d3-7a845f58efdd"
version = "0.108.5+4"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "c8abc88faa3f7a3950832ac5d6e690881590d6dc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[PowerModelsDistribution]]
deps = ["Dates", "InfrastructureModels", "JSON", "JuMP", "LinearAlgebra", "Logging", "LoggingExtras", "MathOptInterface", "Statistics"]
git-tree-sha1 = "2fa6116c0681c76ae6cf8d70527c9fdede9cdc8e"
uuid = "d7431456-977f-11e9-2de3-97ff7677985e"
version = "0.11.4"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "b3fb709f3c97bfc6e948be68beeecb55a0b340ae"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.1"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SCS]]
deps = ["BinaryProvider", "Libdl", "LinearAlgebra", "MathOptInterface", "MathProgBase", "Requires", "SCS_GPU_jll", "SCS_jll", "SparseArrays"]
git-tree-sha1 = "efc1cf417dffb455b74a3d48857ae418098fc8a1"
uuid = "c946c3f1-0d1f-5ce8-9dea-7daa1f7e2d13"
version = "0.7.1"

[[SCS_GPU_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "9c0b47a387c4c7d8a66abb11977391a6602d2384"
uuid = "af6e375f-46ec-5fa0-b791-491b0dfa44a4"
version = "2.1.3+0"

[[SCS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "4816f8272a3e9d08a32a4e94b367c7e84057c145"
uuid = "f4f2fc5b-1d94-523c-97ea-2ab488bedf4b"
version = "2.1.3+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "a50550fa3164a8c46747e62063b4d774ac1bcf49"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.5.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "745914ebcd610da69f3cb6bf76cb7bb83dcb8c9a"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.4"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[Syslogs]]
deps = ["Printf", "Sockets"]
git-tree-sha1 = "46badfcc7c6e74535cc7d833a91f4ac4f805f86d"
uuid = "cea106d9-e007-5e6c-ad93-58fe2094e9c4"
version = "0.3.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimeZones]]
deps = ["Dates", "EzXML", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "960099aed321e05ac649c90d583d59c9309faee1"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.5"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "7c53c35547de1c5b9d46a4797cf6d8253807108c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.5"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[ZipFile]]
deps = ["Libdl", "Printf", "Zlib_jll"]
git-tree-sha1 = "c3a5637e27e914a7a445b8d0ad063d701931e9f7"
uuid = "a5390f91-8eb1-5f08-bee0-b1d1ffed6cea"
version = "0.9.3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─f9e56b1c-a6c4-4704-8d37-151caa33a1cc
# ╠═57d10945-b2b3-4c2f-8ae2-16db9a488ee6
# ╟─00daa957-0a74-4609-9a1b-ea2b788e45a3
# ╟─e8d28cbb-69fb-44de-81f5-86dc903df229
# ╠═ab88f189-3093-4a54-80ef-9120e9a20202
# ╟─79e5480a-aecb-49fb-ba7e-d1904a13b58d
# ╟─c139b036-c4f1-4a1e-aed1-5fc77e2ab00c
# ╟─67033525-d437-4d12-bfe9-c7acb6df1bc2
# ╟─2d0f68f9-c336-4120-85bb-e781b55bcb5f
# ╟─d8749fe0-3915-4b63-a9e3-92726b590962
# ╟─f22095fa-0be0-4b78-86f7-628a015f91f8
# ╟─8f4fd610-d4b7-45bf-9473-65f7720cb224
# ╠═01502fce-b6b9-4696-8562-8f674f1575f8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
