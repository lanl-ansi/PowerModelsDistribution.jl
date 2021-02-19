"multiconductor version of instantiate model from PowerModels"
function instantiate_mc_model(data::Dict{String,<:Any}, model_type::Type, build_method::Function; ref_extensions::Vector{<:Function}=Vector{Function}([]), kwargs...)
    if get(data, "data_model", MATHEMATICAL) == ENGINEERING
        Memento.info(_LOGGER, "Converting ENGINEERING data model to MATHEMATICAL first to build JuMP model")
        data = transform_data_model(data)
    end

    return _IM.instantiate_model(
        data, model_type, build_method, _PM.ref_add_core!, _PM._pm_global_keys,
        pmd_it_sym; ref_extensions = [ref_extensions..., ref_add_arcs_switch!,
        ref_add_arcs_transformer!, ref_add_connections!], kwargs...)
end


# Helper functions for multinetwork AbstractPowerModel objects.
ismultinetwork(pm::_PM.AbstractPowerModel) = _IM.ismultinetwork(pm, pmd_it_sym)
nw_ids(pm::_PM.AbstractPowerModel) = _IM.nw_ids(pm, pmd_it_sym)
nws(pm::_PM.AbstractPowerModel) = _IM.nws(pm, pmd_it_sym)


# Helper functions for AbstractPowerModel component indices.
ids(pm::_PM.AbstractPowerModel, nw::Int, key::Symbol) = _IM.ids(pm, pmd_it_sym, nw, key)
ids(pm::_PM.AbstractPowerModel, key::Symbol; nw::Int=nw_id_default) = _IM.ids(pm, pmd_it_sym, key; nw = nw)


# Helper functions for AbstractPowerModel `ref` access.
ref(pm::_PM.AbstractPowerModel, nw::Int = nw_id_default) = _IM.ref(pm, pmd_it_sym, nw)
ref(pm::_PM.AbstractPowerModel, nw::Int, key::Symbol) = _IM.ref(pm, pmd_it_sym, nw, key)
ref(pm::_PM.AbstractPowerModel, nw::Int, key::Symbol, idx::Int) = _IM.ref(pm, pmd_it_sym, nw, key, idx)
ref(pm::_PM.AbstractPowerModel, nw::Int, key::Symbol, idx::Int, param::String) = _IM.ref(pm, pmd_it_sym, nw, key, idx, param)
ref(pm::_PM.AbstractPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.ref(pm, pmd_it_sym, key; nw = nw)
ref(pm::_PM.AbstractPowerModel, key::Symbol, idx::Int; nw::Int = nw_id_default) = _IM.ref(pm, pmd_it_sym, key, idx; nw = nw)
ref(pm::_PM.AbstractPowerModel, key::Symbol, idx::Int, param::String; nw::Int = nw_id_default) = _IM.ref(pm, pmd_it_sym, key, idx, param; nw = nw)


# Helper functions for AbstractPowerModel `var` access.
var(pm::_PM.AbstractPowerModel, nw::Int = nw_id_default) = _IM.var(pm, pmd_it_sym, nw)
var(pm::_PM.AbstractPowerModel, nw::Int, key::Symbol) = _IM.var(pm, pmd_it_sym, nw, key)
var(pm::_PM.AbstractPowerModel, nw::Int, key::Symbol, idx::Int) = _IM.var(pm, pmd_it_sym, nw, key, idx)
var(pm::_PM.AbstractPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.var(pm, pmd_it_sym, key; nw = nw)
var(pm::_PM.AbstractPowerModel, key::Symbol, idx::Int; nw::Int = nw_id_default) = _IM.var(pm, pmd_it_sym, key, idx; nw = nw)


# Helper functions for AbstractPowerModel `con` access.
con(pm::_PM.AbstractPowerModel, nw::Int = nw_id_default) = _IM.con(pm, pmd_it_sym; nw = nw)
con(pm::_PM.AbstractPowerModel, nw::Int, key::Symbol) = _IM.con(pm, pmd_it_sym, nw, key)
con(pm::_PM.AbstractPowerModel, nw::Int, key::Symbol, idx) = _IM.con(pm, pmd_it_sym, nw, key, idx)
con(pm::_PM.AbstractPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.con(pm, pmd_it_sym, key; nw = nw)
con(pm::_PM.AbstractPowerModel, key::Symbol, idx; nw::Int = nw_id_default) = _IM.con(pm, pmd_it_sym, key, idx; nw = nw)


# Helper functions for AbstractPowerModel `sol` access.
sol(pm::_PM.AbstractPowerModel, nw::Int = nw_id_default) = _IM.sol(pm, pmd_it_sym; nw = nw)
sol(pm::_PM.AbstractPowerModel, nw::Int, key::Symbol) = _IM.sol(pm, pmd_it_sym, nw, key)
sol(pm::_PM.AbstractPowerModel, nw::Int, key::Symbol, idx) = _IM.sol(pm, pmd_it_sym, nw, key, idx)
sol(pm::_PM.AbstractPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.sol(pm, pmd_it_sym, key; nw = nw)
sol(pm::_PM.AbstractPowerModel, key::Symbol, idx; nw::Int = nw_id_default) = _IM.sol(pm, pmd_it_sym, key, idx; nw = nw)