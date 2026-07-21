using PowerModelsDistribution
const PMD = PowerModelsDistribution

import InfrastructureModels

PowerModelsDistribution.silence!()

import JuMP
import Ipopt
import SCS

import JSON

using Test
using LinearAlgebra

using Logging
# global_logger(ConsoleLogger(stderr, Logging.Debug))


pmd_path = joinpath(dirname(pathof(PowerModelsDistribution)), "..")

ipopt_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes")
ipopt_solver_adaptive = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "mu_strategy"=>"adaptive", "sb"=>"yes","warm_start_init_point"=>"yes")
scs_solver = optimizer_with_attributes(SCS.Optimizer, "verbose"=>0)

include("common.jl")
include("test_cases.jl")

@testset "PowerModelsDistribution" begin
    include("opf_ravens.jl")

    include("multinetwork.jl")

    include("opendss.jl")

    include("data.jl")

    include("pf.jl")

    include("pf_bf.jl")

    include("opf.jl")

    include("opf_bf.jl")

    include("opf_iv.jl")

    include("storage.jl")

    include("debug.jl")

    include("transformer.jl")

    include("capacitor.jl")

    include("loadmodels.jl")

    include("delta_gens.jl")

    include("shunt.jl")

    include("mld.jl")

    include("data_model.jl")

    include("en_opf_bounds.jl")

    include("en_pf_validation.jl")

    include("en_pf_native_validation.jl")

    include("line_constants.jl")

    # NEW TESTS

    include("ravens/data.jl")

    include("ravens/pf.jl")

    include("ravens/pf_bf.jl")

    include("ravens/opf.jl")

    include("ravens/opf_bf.jl")

    include("ravens/opf_iv.jl")

    include("ravens/storage.jl")

    include("ravens/debug.jl")

    include("ravens/multinetwork.jl")

    include("ravens/transformer.jl")

    include("ravens/capacitor.jl")

    include("ravens/loadmodels.jl")

    include("ravens/shunt.jl")

    include("ravens/mld.jl")

    include("ravens/opf_ravens_extra.jl")
end
