using PowerModelsDistribution
const PMD = PowerModelsDistribution

import InfrastructureModels

import PowerModels
const PM = PowerModels

PowerModels.silence()
PowerModelsDistribution.silence!()

import JuMP
import Ipopt
import Cbc
import Juniper
import SCS

import JSON

using Test
using LinearAlgebra

pms_path = joinpath(dirname(pathof(PowerModels)), "..")
pmd_path = joinpath(dirname(pathof(PowerModelsDistribution)), "..")

ipopt_solver = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)
ipopt_solver_adaptive = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0, "mu_strategy"=>"adaptive")
cbc_solver = optimizer_with_attributes(Cbc.Optimizer, "logLevel"=>0)
scs_solver = optimizer_with_attributes(SCS.Optimizer, "max_iters"=>20000, "eps"=>1e-5, "alpha"=>0.4, "verbose"=>0)
juniper_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-4, "print_level"=>0), "mip_solver"=>cbc_solver, "log_levels"=>[])

include("common.jl")

@testset "PowerModelsDistribution" begin

    include("opendss.jl")

    include("data.jl")

    include("pf.jl")

    include("pf_bf.jl")

    include("opf.jl")

    include("opf_bf.jl")

    include("opf_iv.jl")

    include("storage.jl")

    include("debug.jl")

    include("multinetwork.jl")

    include("transformer.jl")

    include("loadmodels.jl")

    include("delta_gens.jl")

    include("shunt.jl")

    include("mld.jl")

    include("data_model.jl")

    include("en_opf_bounds.jl")

    include("en_pf_validation.jl")
end
