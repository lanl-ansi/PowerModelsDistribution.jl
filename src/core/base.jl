"root of the power formulation type hierarchy"
abstract type AbstractMCPowerModel <: _IM.AbstractInfrastructureModel end

"a macro for adding the base PowerModels fields to a type definition"
_IM.@def pmd_fields begin
    # this must be explicitly qualified, so that it works in downstream
    # packages that use import PowerModels and this command appears in the
    # downstream package's scope
    PowerModelsDistribution.@im_fields
end


""
function run_model(data::Dict{String,<:Any}, model_type::Type, optimizer, build_method;
        ref_extensions=[], solution_processors=[], relax_integrality=false,
        multinetwork=false, multiconductor=false, kwargs...)

    if multinetwork != _IM.ismultinetwork(data)
        model_requirement = multinetwork ? "multi-network" : "single-network"
        data_type = _IM.ismultinetwork(data) ? "multi-network" : "single-network"
        Memento.error(_LOGGER, "attempted to build a $(model_requirement) model with $(data_type) data")
    end

    if multiconductor != ismulticonductor(data)
        model_requirement = multiconductor ? "multi-conductor" : "single-conductor"
        data_type = ismulticonductor(data) ? "multi-conductor" : "single-conductor"
        Memento.error(_LOGGER, "attempted to build a $(model_requirement) model with $(data_type) data")
    end

    start_time = time()
    pm = instantiate_mc_model(data, model_type, build_method; ref_extensions=ref_extensions, kwargs...)
    Memento.debug(_LOGGER, "pm model build time: $(time() - start_time)")

    start_time = time()
    result = optimize_model!(pm, relax_integrality=relax_integrality, optimizer=optimizer, solution_processors=solution_processors)
    Memento.debug(_LOGGER, "pm model solve and solution time: $(time() - start_time)")

    return result
end


"multiconductor version of instantiate model from PowerModels"
function instantiate_mc_model(data::Dict{String,<:Any}, model_type::Type, build_method::Function; ref_extensions::Vector{<:Function}=Vector{Function}([]), kwargs...)
    if get(data, "data_model", MATHEMATICAL) == ENGINEERING
        Memento.info(_LOGGER, "Converting ENGINEERING data model to MATHEMATICAL first to build JuMP model")
        data = transform_data_model(data)
    end

    return _IM.instantiate_model(
        data, model_type, build_method, ref_add_core!, _pmd_global_keys,
        pmd_it_sym; ref_extensions = [ref_extensions..., ref_add_arcs_switch!,
        ref_add_arcs_transformer!, ref_add_connections!], kwargs...)
end


""
ismulticonductor(pm::AbstractMCPowerModel, nw::Int) = haskey(pm.ref[:it][pmd_it_sym][:nw][nw], :conductors)
ismulticonductor(pm::AbstractMCPowerModel; nw::Int=nw_id_default) = haskey(pm.ref[:it][pmd_it_sym][:nw][nw], :conductors)


# Helper functions for multinetwork AbstractPowerModel objects.
ismultinetwork(pm::AbstractMCPowerModel) = _IM.ismultinetwork(pm, pmd_it_sym)
nw_ids(pm::AbstractMCPowerModel) = _IM.nw_ids(pm, pmd_it_sym)
nws(pm::AbstractMCPowerModel) = _IM.nws(pm, pmd_it_sym)


# Helper functions for AbstractPowerModel component indices.
ids(pm::AbstractMCPowerModel, nw::Int, key::Symbol) = _IM.ids(pm, pmd_it_sym, nw, key)
ids(pm::AbstractMCPowerModel, key::Symbol; nw::Int=nw_id_default) = _IM.ids(pm, pmd_it_sym, key; nw = nw)


# Helper functions for AbstractPowerModel `ref` access.
ref(pm::AbstractMCPowerModel, nw::Int = nw_id_default) = _IM.ref(pm, pmd_it_sym, nw)
ref(pm::AbstractMCPowerModel, nw::Int, key::Symbol) = _IM.ref(pm, pmd_it_sym, nw, key)
ref(pm::AbstractMCPowerModel, nw::Int, key::Symbol, idx::Any) = _IM.ref(pm, pmd_it_sym, nw, key, idx)
ref(pm::AbstractMCPowerModel, nw::Int, key::Symbol, idx::Any, param::String) = _IM.ref(pm, pmd_it_sym, nw, key, idx, param)
ref(pm::AbstractMCPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.ref(pm, pmd_it_sym, key; nw = nw)
ref(pm::AbstractMCPowerModel, key::Symbol, idx::Any; nw::Int = nw_id_default) = _IM.ref(pm, pmd_it_sym, key, idx; nw = nw)
ref(pm::AbstractMCPowerModel, key::Symbol, idx::Any, param::String; nw::Int = nw_id_default) = _IM.ref(pm, pmd_it_sym, key, idx, param; nw = nw)


# Helper functions for AbstractPowerModel `var` access.
var(pm::AbstractMCPowerModel, nw::Int = nw_id_default) = _IM.var(pm, pmd_it_sym, nw)
var(pm::AbstractMCPowerModel, nw::Int, key::Symbol) = _IM.var(pm, pmd_it_sym, nw, key)
var(pm::AbstractMCPowerModel, nw::Int, key::Symbol, idx::Any) = _IM.var(pm, pmd_it_sym, nw, key, idx)
var(pm::AbstractMCPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.var(pm, pmd_it_sym, key; nw = nw)
var(pm::AbstractMCPowerModel, key::Symbol, idx::Any; nw::Int = nw_id_default) = _IM.var(pm, pmd_it_sym, key, idx; nw = nw)


# Helper functions for AbstractPowerModel `con` access.
con(pm::AbstractMCPowerModel, nw::Int = nw_id_default) = _IM.con(pm, pmd_it_sym; nw = nw)
con(pm::AbstractMCPowerModel, nw::Int, key::Symbol) = _IM.con(pm, pmd_it_sym, nw, key)
con(pm::AbstractMCPowerModel, nw::Int, key::Symbol, idx::Any) = _IM.con(pm, pmd_it_sym, nw, key, idx)
con(pm::AbstractMCPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.con(pm, pmd_it_sym, key; nw = nw)
con(pm::AbstractMCPowerModel, key::Symbol, idx::Any; nw::Int = nw_id_default) = _IM.con(pm, pmd_it_sym, key, idx; nw = nw)


# Helper functions for AbstractPowerModel `sol` access.
sol(pm::AbstractMCPowerModel, nw::Int = nw_id_default) = _IM.sol(pm, pmd_it_sym; nw = nw)
sol(pm::AbstractMCPowerModel, nw::Int, key::Symbol) = _IM.sol(pm, pmd_it_sym, nw, key)
sol(pm::AbstractMCPowerModel, nw::Int, key::Symbol, idx::Any) = _IM.sol(pm, pmd_it_sym, nw, key, idx)
sol(pm::AbstractMCPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.sol(pm, pmd_it_sym, key; nw = nw)
sol(pm::AbstractMCPowerModel, key::Symbol, idx::Any; nw::Int = nw_id_default) = _IM.sol(pm, pmd_it_sym, key, idx; nw = nw)
