isdefined(Base, :__precompile__) && __precompile__()

module ThreePhasePowerModels

using JuMP
using PowerModels
PMs = PowerModels

const LOGGER = getlogger(PowerModels)
setlevel!(LOGGER, "info")

end
