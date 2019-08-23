using PowerModelsDistribution
const PMD = PowerModelsDistribution

import Memento

import InfrastructureModels

import PowerModels
const PMs = PowerModels

# Suppress warnings during testing.
const TESTLOG = Memento.getlogger(PowerModels)
Memento.setlevel!(TESTLOG, "error")

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

ipopt_solver = with_optimizer(Ipopt.Optimizer, tol=1e-6, print_level=0)
cbc_solver = with_optimizer(Cbc.Optimizer, logLevel=0)
scs_solver = with_optimizer(SCS.Optimizer, max_iters=10000, verbose=0)
juniper_solver = with_optimizer(Juniper.Optimizer, nl_solver=with_optimizer(Ipopt.Optimizer, tol=1e-4, print_level=0), mip_solver=cbc_solver, log_levels=[])

include("common.jl")

@testset "PowerModelsDistribution" begin

    include("matlab.jl")

    include("opendss.jl")

    include("data.jl")

    include("base.jl")

    include("tp_pf.jl")

    include("tp_opf.jl")

    include("tp_opf_bf.jl")

    include("tp_opf-var.jl")

    include("tp_debug.jl")

    include("tp_multinetwork.jl")

    include("transformer.jl")

    include("loadmodels.jl")

    include("tp_mld.jl")
end
