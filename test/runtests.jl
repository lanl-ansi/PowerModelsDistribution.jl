using ThreePhasePowerModels
using Memento

# Suppress warnings during testing.
setlevel!(getlogger(PowerModels), "error")

using Ipopt

using Base.Test

@testset "ThreePhasePowerModels" begin

include("matlab.jl")

end
