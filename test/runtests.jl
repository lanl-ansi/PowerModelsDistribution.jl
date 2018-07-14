using ThreePhasePowerModels
const TPPMs = ThreePhasePowerModels

using Memento

using PowerModels
const PMs = PowerModels

# Suppress warnings during testing.
setlevel!(getlogger(PowerModels), "error")

using Ipopt
using Cbc
using Pajarito
using Juniper

using Base.Test


ipopt_solver = IpoptSolver(tol=1e-6, print_level=0)
cbc_solver = CbcSolver()
juniper_solver = JuniperSolver(IpoptSolver(tol=1e-4, print_level=0), mip_solver=cbc_solver, log_levels=[])
pajarito_solver = PajaritoSolver(mip_solver=cbc_solver, cont_solver=ipopt_solver, log_level=0)


@testset "ThreePhasePowerModels" begin

    include("matlab.jl")

    include("opendss.jl")

    include("data.jl")

    include("base.jl")

    include("tp_opf.jl")

    # include("tp_ots.jl")

    include("tp_pf.jl")

    include("tp_opf_bf.jl")

    include("tp_debug.jl")

end
