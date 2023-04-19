abstract type MathBus <: MathNodeObject end

"""
"""
Base.@kwdef mutable struct MathBusObj <: MathBus
    id::Int
    bus_type::BusType = PQ
    terminals::Vector{Int}
    gounded::Vector{Bool} = fill(false, length(terminals))
    vm::Vector{Float64}
    va::Vector{Float64}
    vm_lb::Vector{Float64} = fill(0.0, length(terminals))
    vm_ub::Vector{Float64} = fill(Inf, length(terminals))
    vm_pair_lb::Vector{Tuple} = Vector{Tuple}()
    vm_pair_ub::Vector{Tuple} = Vector{Tuple}()
    status::Status = ENABLED
    is_per_unit::Bool = false
    source_id::String
end


abstract type MathBranch <: MathEdgeObject end


"""
"""
Base.@kwdef mutable struct MathBranchObj <: MathBranch
    id::Int
    f_bus::Int
    t_bus::Int
    f_connections::Vector{Int}
    t_connections::Vector{Int}
    r::Matrix{Float64} = fill(0.0, length(f_connections), length(t_connections))
    x::Matrix{Float64} = fill(0.0, length(f_connections), length(t_connections))
    g_fr::Matrix{Float64} = fill(0.0, length(f_connections), length(t_connections))
    b_fr::Matrix{Float64} = fill(0.0, length(f_connections), length(t_connections))
    g_to::Matrix{Float64} = fill(0.0, length(f_connections), length(t_connections))
    b_to::Matrix{Float64} = fill(0.0, length(f_connections), length(t_connections))
    cm_ub::Vector{Float64} = fill(Inf, length(f_connections))
    sm_ub::Vector{Float64} = fill(Inf, length(f_connections))
    vad_lb::Vector{Float64} = fill(-15.0, length(f_connections))
    vad_ub::Vector{Float64} = fill( 15.0, length(f_connections))
    status::Status = ENABLED
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing,DssLine,DssCapacitor,DssReactor} = missing
end


abstract type MathSwitch <: MathEdgeObject end


"""
"""
Base.@kwdef mutable struct MathSwitchObj <: MathSwitch
    id::Int
    f_bus::Int
    t_bus::Int
    f_connections::Vector{Int}
    t_connections::Vector{Int}
    cm_ub::Vector{Float64} = fill(Inf, length(f_connections))
    sm_ub::Vector{Float64} = fill(Inf, length(f_connections))
    state::SwitchState = CLOSED
    dispatchable::Dispatchable = YES
    status::Status = ENABLED
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing,DssLine} = missing
end


abstract type MathTransformer <: MathEdgeObject end

"""
"""
Base.@kwdef mutable struct MathTransformerObj <: MathTransformer
    id::Int
    f_bus::Int
    t_bus::Int
    f_connections::Vector{Int}
    t_connections::Vector{Int}
    configuration::ConnConfig
    tm_nom::Vector{Float64}
    tm_set::Vector{Float64}
    tm_fix::Vector{Bool}
    polarity::Int
    cm_ub::Float64 = Inf
    sm_ub::Float64 = Inf
    status::Status = ENABLED
    controls::Union{EngTransformerControls,Missing} = missing
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing,DssTransformer} = missing
end

abstract type MathLoad <: MathNodeObject end

"""
"""
Base.@kwdef mutable struct MathLoadObj <: MathLoad
    id::Int
    bus::Int
    connections::Vector{Int}
    configuration::ConnConfig = WYE
    model::LoadModel = POWER
    pd::Vector{Float64}
    qd::Vector{Float64}
    vm_nom::Float64
    dispatchable::Dispatchable = NO
    status::Status = ENABLED
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing,DssLoad} = missing
end


abstract type MathShunt <: MathNodeObject end


"""
"""
Base.@kwdef mutable struct MathShuntObj <: MathShunt
    id::Int
    bus::Int
    connections::Vector{Int}
    configuration::ConnConfig = WYE
    gs::Vector{Float64}
    bs::Vector{Float64}
    dispatchable::Dispatchable = NO
    status::Status = ENABLED
    controls::Union{Missing,EngShuntControls} = missing
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing,DssCapacitor} = missing
end

abstract type MathGenerator <: MathNodeObject end

"""
"""
Base.@kwdef mutable struct MathGeneratorObj <: MathGenerator
    id::Int
    bus::Int
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
    dispatchable::Dispatchable = YES
    status::Status = ENABLED
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing, DssGenerator, DssPvsystem, DssVsource}
end


abstract type MathStorage <: MathNodeObject end

"""
"""
Base.@kwdef mutable struct MathStorageObj <: MathStorage
    id::Int
    bus::Int
    connections::Vector{Int}
    control_mode::ControlMode = FREQUENCYDROOP
    configuration::ConnConfig = WYE
    energy::Float64
    energy_ub::Float64
    charge_ub::Float64
    discharge_ub::Float64
    charge_efficiency::Float64
    discharge_efficiency::Float64
    sm_ub::Union{Vector{Float64},Float64}
    qs_lb::Union{Vector{Float64},Float64}
    qs_ub::Union{Vector{Float64},Float64}
    r::Union{Vector{Float64},Float64}
    x::Union{Vector{Float64},Float64}
    p_loss::Union{Vector{Float64},Float64}
    q_loss::Union{Vector{Float64},Float64}
    ps::Union{Vector{Float64},Float64}
    qs::Union{Vector{Float64},Float64}
    dispatchable::Dispatchable = YES
    status::Status = ENABLED
    is_per_unit::Bool = false
    source_id::String
    dss::Union{Missing,DssStorage}
end


"""
"""
Base.@kwdef struct MathematicalDataModel <: MathematicalModel{SubnetworkModel}
    settings::Settings = Settings()
    metadata::Metadata = Metadata()
    bus_lookup::Dict{String,Int} = Dict{String,Int}()

    bus::Dict{Int,<:MathBus} = Dict{Int,MathBus}()

    branch::Dict{Int,<:MathBranch} = Dict{Int,MathBranch}()
    switch::Dict{Int,<:MathSwitch} = Dict{Int,MathSwitch}()
    transformer::Dict{Int,<:MathTransformer} = Dict{Int,MathTransformer}()

    load::Dict{Int,<:MathLoad} = Dict{Int,MathLoad}()
    shunt::Dict{Int,<:MathShunt} = Dict{Int,MathShunt}()
    generator::Dict{Int,<:MathGenerator} = Dict{Int,MathGenerator}()
    storage::Dict{Int,<:MathStorage} = Dict{Int,MathStorage}()
end


"""
"""
Base.@kwdef struct MultinetworkMathematicalDataModel <: MathematicalModel{MultinetworkModel}
    metadata::Metadata = Metadata()

    nw::Dict{Int,<:MathematicalModel} = Dict{Int,MathematicalModel}()
    nw_map::Dict{String,Int}
end
