module ThreePhasePowerModels

import JuMP
import PowerModels
import InfrastructureModels
import Memento

import LinearAlgebra

const PMs = PowerModels

function __init__()
    global LOGGER = Memento.getlogger(PowerModels)
    Memento.warn(LOGGER, "DEPRECIATED :: ThreePhasePowerModels.jl is depreciated, v0.3.2 is the final version with this package name. PowerModelsDistribution.jl (https://github.com/lanl-ansi/PowerModelsDistribution.jl) is the new name of this package going forward.")
end

include("core/ref.jl")
include("core/multiconductor.jl")
include("core/variable.jl")
include("core/constraint.jl")
include("core/objective.jl")
include("core/solution.jl")

include("form/acp.jl")
include("form/acr.jl")
include("form/apo.jl")
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
include("prob/tp_opf_lm.jl")
include("prob/tp_opf_oltc.jl")
include("prob/tp_opf_bctr.jl")
include("prob/tp_opf_bf.jl")
include("prob/tp_ots.jl")
include("prob/tp_pf.jl")
include("prob/tp_pf_lm.jl")
include("prob/tp_pf_bf.jl")
include("prob/tp_debug.jl")
include("prob/tp_test.jl")

end
