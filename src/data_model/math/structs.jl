abstract type MathBus <: MathNodeObject end

"""
"""
Base.@kwdef mutable struct MathBusObj <: MathBus
    index::Int
    bus_i::Int
    bus_type::Int = 1
    terminals::Vector{Int}
    gounded::Vector{Bool} = fill(false, length(terminals))
    vm::Vector{Float64}
    va::Vector{Float64}
    vmin::Vector{Float64} = fill(0.0, length(terminals))
    vmax::Vector{Float64} = fill(Inf, length(terminals))
    vm_pair_lb::Vector{Tuple} = Vector{Tuple}()
    vm_pair_ub::Vector{Tuple} = Vector{Tuple}()
    vm_pn_lb::Union{Missing,Float64} = missing
    vm_pn_ub::Union{Missing,Float64} = missing
    vm_pp_lb::Union{Missing,Float64} = missing
    vm_pp_ub::Union{Missing,Float64} = missing
    vm_ng_lb::Union{Missing,Float64} = missing
    vm_ng_ub::Union{Missing,Float64} = missing
    is_per_unit::Bool = false
    source_id::String
end


abstract type MathBranch <: MathEdgeObject end


"""
"""
Base.@kwdef mutable struct MathBranchObj <: MathBranch
    index::Int
    f_bus::Int
    t_bus::Int
    f_connections::Vector{Int}
    t_connections::Vector{Int}
    br_r::Matrix{Float64} = fill(0.0, length(f_connections), length(t_connections))
    br_x::Matrix{Float64} = fill(0.0, length(f_connections), length(t_connections))
    g_fr::Matrix{Float64} = fill(0.0, length(f_connections), length(t_connections))
    b_fr::Matrix{Float64} = fill(0.0, length(f_connections), length(t_connections))
    g_to::Matrix{Float64} = fill(0.0, length(f_connections), length(t_connections))
    b_to::Matrix{Float64} = fill(0.0, length(f_connections), length(t_connections))
    c_rating_a::Vector{Float64} = fill(Inf, length(f_connections))
    c_rating_b::Vector{Float64} = fill(Inf, length(f_connections))
    c_rating_c::Vector{Float64} = fill(Inf, length(f_connections))
    rate_a::Vector{Float64} = fill(Inf, length(f_connections))
    rate_b::Vector{Float64} = fill(Inf, length(f_connections))
    rate_c::Vector{Float64} = fill(Inf, length(f_connections))
    angmin::Vector{Float64} = fill(-15.0, length(f_connections))
    angmax::Vector{Float64} = fill( 15.0, length(f_connections))
    br_status::Int = 1
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing,DssLine,DssCapacitor,DssReactor} = missing
end


abstract type MathSwitch <: MathEdgeObject end


"""
"""
Base.@kwdef mutable struct MathSwitchObj <: MathSwitch
    index::Int
    f_bus::Int
    t_bus::Int
    f_connections::Vector{Int}
    t_connections::Vector{Int}
    current_rating::Vector{Float64} = fill(Inf, length(f_connections))
    c_rating_b::Vector{Float64} = fill(Inf, length(f_connections))
    c_rating_c::Vector{Float64} = fill(Inf, length(f_connections))
    thermal_rating::Vector{Float64} = fill(Inf, length(f_connections))
    rate_b::Vector{Float64} = fill(Inf, length(f_connections))
    rate_c::Vector{Float64} = fill(Inf, length(f_connections))
    state::Int = 1
    dispatchable::Int = 0
    status::Int = 1
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing,DssLine} = missing
end

abstract type MathTransformerControls <: MathDataObject end

Base.@kwdef mutable struct MathTransformerControlsObj <: MathTransformerControls
    vreg
    band
    ptratio
    ctprim
    r
    x
end

abstract type MathTransformer <: MathEdgeObject end

"""
"""
Base.@kwdef mutable struct MathTransformerObj <: MathTransformer
    index::Int
    f_bus::Int
    t_bus::Int
    f_connections::Vector{Int}
    t_connections::Vector{Int}
    configuration::ConnConfig
    tm_nom::Float64 = 1.0
    tm_set::Vector{Float64} = fill(1.0, length(f_connections))
    tm_fix::Vector{Bool} = fill(true, length(f_connections))
    polarity::Int = -1
    cm_ub::Float64 = Inf
    sm_ub::Float64 = Inf
    status::Int = 1
    controls::Union{MathTransformerControls,Missing} = missing
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing,DssTransformer} = missing
end

abstract type MathLoad <: MathNodeObject end

"""
"""
Base.@kwdef mutable struct MathLoadObj <: MathLoad
    index::Int
    load_bus::Int
    connections::Vector{Int}
    configuration::ConnConfig = WYE
    model::LoadModel = POWER
    pd::Vector{Float64}
    qd::Vector{Float64}
    vnom_kv::Float64
    dispatchable::Dispatchable = NO
    status::Status = ENABLED
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing,DssLoad} = missing
end


abstract type MathShuntControls <: MathDataObject end


"""
"""
Base.@kwdef mutable struct MathShuntControlsObj <: MathShuntControls
    element::Dict{String,Any} = Dict{String,Any}()
end


abstract type MathShunt <: MathNodeObject end


"""
"""
Base.@kwdef mutable struct MathShuntObj <: MathShunt
    index::Int
    shunt_bus::Int
    connections::Vector{Int}
    configuration::ConnConfig = WYE
    gs::Vector{Float64}
    bs::Vector{Float64}
    dispatchable::Int = 1
    status::Int = 1
    controls::Union{Missing,MathShuntControlsObj} = missing
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing,DssCapacitor} = missing
end

abstract type MathGenerator <: MathNodeObject end

"""
"""
Base.@kwdef mutable struct MathGeneratorObj <: MathGenerator
    index::Int
    gen_bus::Int
    connections::Vector{Int}
    control_mode::ControlMode = FREQUENCYDROOP
    configuration::ConnConfig = WYE
    pg::Union{Vector{Float64},Float64}
    qg::Union{Vector{Float64},Float64}
    vg::Union{Vector{Float64},Float64}
    pg_lb::Union{Vector{Float64},Float64}
    pg_ub::Union{Vector{Float64},Float64}
    qg_lb::Union{Vector{Float64},Float64}
    qg_ub::Union{Vector{Float64},Float64}
    sm_ub::Union{Vector{Float64},Float64}
    cm_ub::Union{Vector{Float64},Float64}
    dispatchable::Int = 1
    gen_status::Int = 1
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing, DssGenerator, DssPvsystem, DssVsource}
end


abstract type MathStorage <: MathNodeObject end

"""
"""
Base.@kwdef mutable struct MathStorageObj <: MathStorage
    index::Int
    storage_bus::Int
    connections::Vector{Int}
    control_mode::ControlMode = FREQUENCYDROOP
    configuration::ConnConfig = WYE
    energy::Float64
    energy_rating::Float64
    charge_rating::Float64
    discharge_rating::Float64
    charge_efficiency::Float64
    discharge_efficiency::Float64
    thermal_rating::Union{Vector{Float64},Float64}
    qmin::Union{Vector{Float64},Float64}
    qmax::Union{Vector{Float64},Float64}
    r::Union{Vector{Float64},Float64}
    x::Union{Vector{Float64},Float64}
    p_loss::Union{Vector{Float64},Float64}
    q_loss::Union{Vector{Float64},Float64}
    ps::Union{Vector{Float64},Float64}
    qs::Union{Vector{Float64},Float64}
    dispatchable::Int = 1
    status::Int = 1
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing,DssStorage}
end


"""
"""
Base.@kwdef struct MathematicalDataModel <: MathematicalModel{NetworkModel}
    settings::Settings = SettingsObj()
    metadata::Metadata = MetadataObj()
    bus_lookup::Dict{String,Int} = Dict{String,Int}()
    # map::Vector{Dict{String,String}} = Dict{String,String}[]

    bus::Dict{String,<:MathBus} = Dict{String,MathBus}()

    branch::Dict{String,<:MathBranch} = Dict{String,MathBranch}()
    switch::Dict{String,<:MathSwitch} = Dict{String,MathSwitch}()
    transformer::Dict{String,<:MathTransformer} = Dict{String,MathTransformer}()

    load::Dict{String,<:MathLoad} = Dict{String,MathLoad}()
    shunt::Dict{String,<:MathShunt} = Dict{String,MathShunt}()
    generator::Dict{String,<:MathGenerator} = Dict{String,MathGenerator}()
    storage::Dict{String,<:MathStorage} = Dict{String,MathStorage}()
end


"""
"""
Base.@kwdef struct MultinetworkMathematicalDataModel <: MathematicalModel{MultinetworkModel}
    metadata::Metadata = MetadataObj()

    nw::Dict{Int,<:MathematicalModel} = Dict{Int,MathematicalModel}()
    nw_map::Dict{String,Int}
end
