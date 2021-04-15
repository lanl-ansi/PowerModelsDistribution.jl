### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ 30417ccf-d0b4-43d5-9ceb-5b00ef9c8a17
begin
	using PowerModelsDistribution
	using LinearAlgebra: diagm
	using Ipopt
end

# ╔═╡ ff4b2ffa-9d25-11eb-07d0-29751e7b57c2
md"""
# Building Engineering Data Models with Helper Functions

In this notebook we will demonstrate an easy way to start building new distribution network models in the engineering format using new helper functions added in PowerModelsDistribution v0.9
"""

# ╔═╡ f2654d95-0c8b-4524-ad1f-f0f12a35c570
md"""
First, we need a optimizer. In this case we will use Ipopt and initialize it with JuMP's `optimizer_with_attributes`, which we have exported from PowerModelsDistribution by default for the user
"""

# ╔═╡ 1885ce0e-fe1a-4b9b-b652-72cbee21339d
ipopt_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "tol"=>1e-6)

# ╔═╡ 7683dabc-772f-4fdf-a784-549c05084958
md"Initialize an empty `ENGINEERING` model"

# ╔═╡ 20d3d95d-7817-45de-89e9-ffce0202d8cf
Model()

# ╔═╡ d4ce49bd-f38e-44e0-aa3a-4ee3f0f70bf1
md"""
In this block we build a three bus network, with neutrals grounded at the source and loads.

We start with buses, with the sourcebus and loadbus having 4 terminals, with the last terminal grounded.

Then we add a generation source, in this case a voltage source, which is `WYE` configured by default, and therefore expects the last conductor to be a grounded neutral

We add two three phase lines to connect the buses `sourcebus`, `primary`, and `loadbus`. Note that none of the lines have a neutral conductor.

We finally add a three-phase load a the `loadbus` bus, but note again that this is a `WYE` configured load, and like the voltage source, this implies that the last conductor is a grounded neutral for the purposes of kron reduction (which is required by default until explicit 4-wire modeling is added to PowerModelsDistribution)
    
Lastly, we need to define the default vbase of the system at the `sourcebus`
"""

# ╔═╡ 1ca9657b-85b4-4ea6-901c-d36e23e8440f
begin
	eng = Model()

	add_bus!(eng, "sourcebus"; terminals=[1,2,3,4], grounded=[4])
	add_bus!(eng, "primary"; terminals=[1,2,3])
	add_bus!(eng, "loadbus"; terminals=[1,2,3,4], grounded=[4])

	add_voltage_source!(eng, "source", "sourcebus", [1,2,3,4]; vm=[1, 1, 1, 0])

	add_linecode!(eng, "default", diagm(0=>fill(0.01, 3)), diagm(0=>fill(0.2, 3)))

	add_line!(eng, "trunk", "sourcebus", "primary", [1,2,3], [1,2,3]; linecode="default")
	add_line!(eng, "primary", "primary", "loadbus", [1,2,3], [1,2,3]; linecode="default")

	add_load!(eng, "balanced", "loadbus", [1,2,3,4]; pd_nom=[5, 5, 5], qd_nom=[1, 1, 1])

	add_vbase_default!(eng, "sourcebus", 1)
	
	eng
end

# ╔═╡ 9010cccb-e174-4dff-b934-7e7be5a93a1d
md"Running this case with OPF gives the results below"

# ╔═╡ 3ed441b8-2055-41b7-a4b5-29130d5880ed
result = solve_mc_opf(eng, ACPUPowerModel, ipopt_solver)

# ╔═╡ 650b6dee-e753-418c-af1e-a0ddda09df81
md"""In the following example, we provide examples of a wider range of component types that can be added using helper functions to give a flavor of what is possible in PowerModelsDistribution v0.9"""

# ╔═╡ 95ecd659-45d0-4bd0-9f32-9cafd78fb2c7
begin
	eng2 = deepcopy(eng)

	add_bus!(eng2, "ttbus"; terminals=[1,2,3,4], grounded=[4])

	add_transformer!(eng2, "tx1", "sourcebus", "ttbus", [1,2,3,4], [1,2,3,4])

	add_bus!(eng2, "loadbus2"; terminals=[1,2,3,4], grounded=[4])

	add_switch!(eng2, "breaker", "ttbus", "loadbus2", [1,2,3], [1,2,3])

	add_load!(eng2, "tload", "loadbus2", [1,2,3,4]; pd_nom=[5, 5, 5], qd_nom=[1, 1, 1])

	add_generator!(eng2, "secondary", "loadbus2", [1,2,3,4]; cost_pg_parameters=[0.0, 1.2, 0])

	add_shunt!(eng2, "cap", "loadbus2", [1,2,3,4]; bs=diagm(0=>fill(1, 3)))


	eng2
end

# ╔═╡ fdaed782-529b-4267-bf73-e26a6a3d4297
result2 = solve_mc_opf(eng2, ACPUPowerModel, ipopt_solver)

# ╔═╡ Cell order:
# ╟─ff4b2ffa-9d25-11eb-07d0-29751e7b57c2
# ╠═30417ccf-d0b4-43d5-9ceb-5b00ef9c8a17
# ╠═f2654d95-0c8b-4524-ad1f-f0f12a35c570
# ╠═1885ce0e-fe1a-4b9b-b652-72cbee21339d
# ╟─7683dabc-772f-4fdf-a784-549c05084958
# ╠═20d3d95d-7817-45de-89e9-ffce0202d8cf
# ╟─d4ce49bd-f38e-44e0-aa3a-4ee3f0f70bf1
# ╠═1ca9657b-85b4-4ea6-901c-d36e23e8440f
# ╟─9010cccb-e174-4dff-b934-7e7be5a93a1d
# ╠═3ed441b8-2055-41b7-a4b5-29130d5880ed
# ╟─650b6dee-e753-418c-af1e-a0ddda09df81
# ╠═95ecd659-45d0-4bd0-9f32-9cafd78fb2c7
# ╠═fdaed782-529b-4267-bf73-e26a6a3d4297
