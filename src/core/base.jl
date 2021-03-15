"root of the power formulation type hierarchy"
abstract type AbstractMCPowerModel <: _IM.AbstractInfrastructureModel end

"a macro for adding the base PowerModels fields to a type definition"
_IM.@def pmd_fields begin
    # this must be explicitly qualified, so that it works in downstream
    # packages that use import PowerModels and this command appears in the
    # downstream package's scope
    PowerModelsDistribution.@im_fields
end


# Helper functions for multinetwork AbstractPowerModel objects.
ismultinetwork(pm::AbstractMCPowerModel) = ismultinetwork(pm, pmd_it_sym)
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
