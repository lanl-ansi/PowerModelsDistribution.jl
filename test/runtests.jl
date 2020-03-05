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
scs_solver = with_optimizer(SCS.Optimizer, max_iters=20000, eps=1e-5, alpha=0.4, verbose=0)
juniper_solver = with_optimizer(Juniper.Optimizer, nl_solver=with_optimizer(Ipopt.Optimizer, tol=1e-4, print_level=0), mip_solver=cbc_solver, log_levels=[])

include("common.jl")

@testset "PowerModelsDistribution" begin

    # include("opendss.jl")

    include("data.jl") # all passing

    include("pf.jl") # all passing

    include("pf_iv.jl") # all passing

    include("opf.jl") # all passing

    include("opf_bf.jl") # all passing

    include("opf_iv.jl") # all passing

    include("storage.jl") # all passing

    include("debug.jl") # all passing

    include("multinetwork.jl") # all passing

    # include("transformer.jl")

    # include("loadmodels.jl")

    # include("delta_gens.jl")

    include("shunt.jl")

    include("mld.jl") # only transformer tests failing
end
