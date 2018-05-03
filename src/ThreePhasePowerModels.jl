# isdefined(Base, :__precompile__) && __precompile__()

module ThreePhasePowerModels

using JuMP
using PowerModels
using InfrastructureModels
using Memento

import Compat: @__MODULE__

const PMs = PowerModels

const LOGGER = getlogger(@__MODULE__)
setlevel!(LOGGER, "info")

include("core/constraint_template.jl")
include("core/ref.jl")
include("core/variable.jl")

include("form/acp.jl")
include("form/dcp.jl")
include("form/df.jl")
include("form/shared.jl")
include("form/wr.jl")

include("io/matlab.jl")
include("io/common.jl")
include("io/opendss.jl")

include("prob/tp_opf.jl")
include("prob/tp_opf_bf.jl")
include("prob/tp_ots.jl")
include("prob/tp_pf.jl")

end
