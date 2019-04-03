module ThreePhasePowerModels

using JuMP
using PowerModels
using InfrastructureModels
using Memento

using Compat.LinearAlgebra

if VERSION < v"0.7.0-"
    import Compat: occursin
    import Compat: findall
    import Compat: undef
    import Compat: Nothing
end

const PMs = PowerModels

function __init__()
    global LOGGER = getlogger(PowerModels)
end

include("core/ref.jl")
include("core/multiconductor.jl")
include("core/variable.jl")
include("core/constraint.jl")
include("core/objective.jl")

include("form/acp.jl")
include("form/acr.jl")
include("form/dcp.jl")
include("form/bf.jl")
include("form/bf_mx.jl")
include("form/bf_mx_lin.jl")
include("form/bf_mx_soc.jl")
include("form/bf_mx_sdp.jl")
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
include("prob/tp_pf_bf.jl")
include("prob/tp_debug.jl")
include("prob/tp_test.jl")

end
