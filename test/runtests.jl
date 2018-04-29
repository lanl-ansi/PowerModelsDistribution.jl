using ThreePhasePowerModels
using Memento

using PowerModels
const PMs = PowerModels

# Suppress warnings during testing.
setlevel!(getlogger(PowerModels), "error")

using Ipopt

using Base.Test


ipopt_solver = IpoptSolver(tol=1e-6, print_level=0)


@testset "ThreePhasePowerModels" begin

include("matlab.jl")

include("data.jl")

include("base.jl")

include("tp_opf.jl")

end
