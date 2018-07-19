# isdefined(Base, :__precompile__) && __precompile__()

module ThreePhasePowerModels

using JuMP
using PowerModels
using InfrastructureModels
using Memento

const PMs = PowerModels

const LOGGER = getlogger(PowerModels)

include("core/ref.jl")
include("core/variable.jl")
include("core/multiphase.jl")
include("core/objective.jl")

include("form/acp.jl")
include("form/dcp.jl")
include("form/bf.jl")
include("form/bf_mat.jl")
include("form/bf_mat_lin.jl")
include("form/bf_mat_soc.jl")
include("form/bf_mat_sdp.jl")
include("form/shared.jl")
include("form/wr.jl")

include("core/constraint_template.jl")
include("core/relaxation_scheme.jl")


include("io/matlab.jl")
include("io/common.jl")
include("io/dss_parse.jl")
include("io/dss_structs.jl")
include("io/opendss.jl")

include("prob/tp_opf.jl")
include("prob/tp_opf_bf.jl")
include("prob/tp_ots.jl")
include("prob/tp_pf.jl")
include("prob/tp_debug.jl")

end
