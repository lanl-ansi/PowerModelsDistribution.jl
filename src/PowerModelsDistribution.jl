module PowerModelsDistribution

    import JSON
    import CSV

    import JuMP
    import MathOptInterface

    import InfrastructureModels

    import Logging
    import LoggingExtras

    import Dates
    import LinearAlgebra

    import LinearAlgebra: diagm
    import Statistics: mean, std

    const _IM = InfrastructureModels

    import InfrastructureModels: optimize_model!, @im_fields, nw_id_default, ismultinetwork, update_data!

    const _pmd_global_keys = Set(["time_series", "per_unit"])
    const pmd_it_name = "pmd"
    const pmd_it_sym = Symbol(pmd_it_name)

    include("core/logging.jl")

    function __init__()
        global _DEFAULT_LOGGER = Logging.current_logger()
        global _LOGGER = Logging.ConsoleLogger(; meta_formatter=PowerModelsDistribution._pmd_metafmt)

        Logging.global_logger(_LOGGER)
    end

    include("core/base.jl")
    include("core/types.jl")
    include("core/data.jl")
    include("core/ref.jl")
    include("core/variable.jl")
    include("core/variable_mx.jl")
    include("core/constraint.jl")
    include("core/objective.jl")
    include("core/solution.jl")

    include("form/acp.jl")
    include("form/acr.jl")
    include("form/apo.jl")
    include("form/dcp.jl")
    include("form/ivr.jl")
    include("form/bf.jl")
    include("form/bf_mx.jl")
    include("form/bf_mx_lin.jl")
    include("form/bf_fbs.jl")
    include("form/fotp.jl")
    include("form/fotr.jl")
    include("form/bf_mx_soc.jl")
    include("form/bf_mx_sdp.jl")
    include("form/shared.jl")
    include("form/wr.jl")

    include("core/constraint_template.jl")
    include("core/relaxation_scheme.jl")

    include("io/utils.jl")
    include("io/dss/dss_parse.jl")
    include("io/dss/dss_data_structs.jl")
    include("io/dss/dss_node_structs.jl")
    include("io/dss/dss_edge_structs.jl")
    include("io/dss/dss_structs.jl")
    include("io/dss/dss2eng.jl")
    include("io/json/json.jl")
    include("io/common.jl")

    include("data_model/utils.jl")
    include("data_model/checks.jl")
    include("data_model/components.jl")
    include("data_model/eng2math.jl")
    include("data_model/math2eng.jl")
    include("data_model/multinetwork.jl")
    include("data_model/transformations.jl")
    include("data_model/units.jl")

    include("prob/common.jl")
    include("prob/mld.jl")
    include("prob/opf.jl")
    include("prob/opf_oltc.jl")
    include("prob/opf_capc.jl")
    include("prob/pf.jl")
    include("prob/debug.jl")
    include("prob/test.jl")
    include("prob/osw.jl")

    include("core/export.jl")

end
