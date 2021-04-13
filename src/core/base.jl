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


"checks if a sufficient number of variables exist for the given keys collection"
function _check_var_keys(vars, keys, var_name, comp_name)
    if length(vars) < length(keys)
        error("$(var_name) decision variables appear to be missing for $(comp_name) components")
    end
end


"detection of whether a constraint should be NL or not"
macro smart_constraint(model, vars, expr)
    esc(quote
        if _has_nl_expression($vars)
            JuMP.@NLconstraint($model, $expr)
        else
            JuMP.@constraint($model, $expr)
        end
    end)
end


"Local wrapper method for JuMP.set_lower_bound, which skips NaN and infinite (-Inf only)"
function set_lower_bound(x::JuMP.VariableRef, bound; loose_bounds::Bool=false, pm=missing, category::Symbol=:default)
    if !(isnan(bound) || bound==-Inf)
        JuMP.set_lower_bound(x, bound)
    elseif loose_bounds
        lbs = pm.ext[:loose_bounds]
        JuMP.set_lower_bound(x, -lbs.bound_values[category])
        push!(lbs.loose_lb_vars, x)
    end
end


"Local wrapper method for JuMP.set_upper_bound, which skips NaN and infinite (+Inf only)"
function set_upper_bound(x::JuMP.VariableRef, bound; loose_bounds::Bool=false, pm=missing, category::Symbol=:default)
    if !(isnan(bound) || bound==Inf)
        JuMP.set_upper_bound(x, bound)
    elseif loose_bounds
        lbs = pm.ext[:loose_bounds]
        JuMP.set_upper_bound(x, lbs.bound_values[category])
        push!(lbs.loose_ub_vars, x)
    end
end


""
function comp_start_value(comp::Dict{String,<:Any}, key::String, conductor::Int, default)
    cond_ind = _get_conductor_indicator(comp)
    if haskey(comp, key) && !isempty(cond_ind)
        return comp[key][findfirst(isequal(conductor), comp[cond_ind])]
    else
        return default
    end
end


""
function comp_start_value(comp::Dict{String,<:Any}, key::String, default=0.0)
    return get(comp, key, default)
end
