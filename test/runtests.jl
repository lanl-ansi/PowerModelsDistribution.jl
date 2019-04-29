using ThreePhasePowerModels
const TPPMs = ThreePhasePowerModels

import Memento

import InfrastructureModels

import PowerModels
const PMs = PowerModels

using JuMP

# Suppress warnings during testing.
const TESTLOG = Memento.getlogger(PowerModels)
Memento.setlevel!(TESTLOG, "error")

import Ipopt
import Cbc
import Pavito
import Juniper
import SCS

using Test
using LinearAlgebra

pms_path = joinpath(dirname(pathof(PowerModels)), "..")

ipopt_solver = Ipopt.IpoptSolver(tol=1e-6, print_level=0)
cbc_solver = Cbc.CbcSolver()
scs_solver = SCS.SCSSolver(max_iters=10000, verbose=0)
juniper_solver = Juniper.JuniperSolver(Ipopt.IpoptSolver(tol=1e-4, print_level=0), mip_solver=cbc_solver, log_levels=[])
pavito_solver = Pavito.PavitoSolver(mip_solver=cbc_solver, cont_solver=ipopt_solver, mip_solver_drives = false, log_level=0)


@testset "TPPMs" begin

    include("matlab.jl")

    include("opendss.jl")

    include("data.jl")

    include("base.jl")

    include("tp_pf.jl")

    include("tp_opf.jl")

    include("tp_opf_bf.jl")

    include("tp_opf-var.jl")

    include("tp_debug.jl")

    ## include("tp_ots.jl")

    include("tp_multinetwork.jl")

    include("transformer.jl")
end
