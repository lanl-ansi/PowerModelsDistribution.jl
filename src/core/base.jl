"root of the power formulation type hierarchy"
abstract type AbstractUnbalancedPowerModel <: _IM.AbstractInfrastructureModel end

"a macro for adding the base PowerModels fields to a type definition"
_IM.@def pmd_fields begin
    # this must be explicitly qualified, so that it works in downstream
    # packages that use import PowerModels and this command appears in the
    # downstream package's scope
    PowerModelsDistribution.@im_fields
end


# Helper functions for multinetwork AbstractPowerModel objects.
nw_ids(pm::AbstractUnbalancedPowerModel) = _IM.nw_ids(pm, pmd_it_sym)
nws(pm::AbstractUnbalancedPowerModel) = _IM.nws(pm, pmd_it_sym)


# Helper functions for AbstractPowerModel component indices.
ids(pm::AbstractUnbalancedPowerModel, nw::Int, key::Symbol) = _IM.ids(pm, pmd_it_sym, nw, key)
ids(pm::AbstractUnbalancedPowerModel, key::Symbol; nw::Int=nw_id_default) = _IM.ids(pm, pmd_it_sym, key; nw = nw)


# Helper functions for AbstractPowerModel `ref` access.
ref(pm::AbstractUnbalancedPowerModel, nw::Int = nw_id_default) = _IM.ref(pm, pmd_it_sym, nw)
ref(pm::AbstractUnbalancedPowerModel, nw::Int, key::Symbol) = _IM.ref(pm, pmd_it_sym, nw, key)
ref(pm::AbstractUnbalancedPowerModel, nw::Int, key::Symbol, idx::Any) = _IM.ref(pm, pmd_it_sym, nw, key, idx)
ref(pm::AbstractUnbalancedPowerModel, nw::Int, key::Symbol, idx::Any, param::String) = _IM.ref(pm, pmd_it_sym, nw, key, idx, param)
ref(pm::AbstractUnbalancedPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.ref(pm, pmd_it_sym, key; nw = nw)
ref(pm::AbstractUnbalancedPowerModel, key::Symbol, idx::Any; nw::Int = nw_id_default) = _IM.ref(pm, pmd_it_sym, key, idx; nw = nw)
ref(pm::AbstractUnbalancedPowerModel, key::Symbol, idx::Any, param::String; nw::Int = nw_id_default) = _IM.ref(pm, pmd_it_sym, key, idx, param; nw = nw)


# Helper functions for AbstractPowerModel `var` access.
var(pm::AbstractUnbalancedPowerModel, nw::Int = nw_id_default) = _IM.var(pm, pmd_it_sym, nw)
var(pm::AbstractUnbalancedPowerModel, nw::Int, key::Symbol) = _IM.var(pm, pmd_it_sym, nw, key)
var(pm::AbstractUnbalancedPowerModel, nw::Int, key::Symbol, idx::Any) = _IM.var(pm, pmd_it_sym, nw, key, idx)
var(pm::AbstractUnbalancedPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.var(pm, pmd_it_sym, key; nw = nw)
var(pm::AbstractUnbalancedPowerModel, key::Symbol, idx::Any; nw::Int = nw_id_default) = _IM.var(pm, pmd_it_sym, key, idx; nw = nw)


# Helper functions for AbstractPowerModel `con` access.
con(pm::AbstractUnbalancedPowerModel, nw::Int = nw_id_default) = _IM.con(pm, pmd_it_sym; nw = nw)
con(pm::AbstractUnbalancedPowerModel, nw::Int, key::Symbol) = _IM.con(pm, pmd_it_sym, nw, key)
con(pm::AbstractUnbalancedPowerModel, nw::Int, key::Symbol, idx::Any) = _IM.con(pm, pmd_it_sym, nw, key, idx)
con(pm::AbstractUnbalancedPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.con(pm, pmd_it_sym, key; nw = nw)
con(pm::AbstractUnbalancedPowerModel, key::Symbol, idx::Any; nw::Int = nw_id_default) = _IM.con(pm, pmd_it_sym, key, idx; nw = nw)


# Helper functions for AbstractPowerModel `sol` access.
sol(pm::AbstractUnbalancedPowerModel, nw::Int = nw_id_default) = _IM.sol(pm, pmd_it_sym; nw = nw)
sol(pm::AbstractUnbalancedPowerModel, nw::Int, key::Symbol) = _IM.sol(pm, pmd_it_sym, nw, key)
sol(pm::AbstractUnbalancedPowerModel, nw::Int, key::Symbol, idx::Any) = _IM.sol(pm, pmd_it_sym, nw, key, idx)
sol(pm::AbstractUnbalancedPowerModel, key::Symbol; nw::Int = nw_id_default) = _IM.sol(pm, pmd_it_sym, key; nw = nw)
sol(pm::AbstractUnbalancedPowerModel, key::Symbol, idx::Any; nw::Int = nw_id_default) = _IM.sol(pm, pmd_it_sym, key, idx; nw = nw)

@doc "helper function to access the `ids` of multinetworks from AbstractUnbalancedPowerModel structs, returns ints" nw_ids
@doc "helper function to access multinetwork data from AbstractUnbalancedPowerModel structs, returns (id,data) pairs" nws
@doc "helper function to access the `ids` of AbstractUnbalancedPowerModel structs' `ref`, returns ints" ids
@doc "helper function to access the AbstractUnbalancedPowerModel structs' `ref`, returns (id,data) pairs" ref
@doc "helper function to access the AbstractUnbalancedPowerModel structs' `var`, returns JuMP VariableRef" var
@doc "helper function to access the AbstractUnbalancedPowerModel structs' `con`, returns JuMP Constraint" con
@doc "helper function to access the AbstractUnbalancedPowerModel structs' `sol`, returns Dict" sol

"checks if a sufficient number of variables exist for the given keys collection"
function _check_var_keys(vars, keys, var_name, comp_name)
    if length(vars) < length(keys)
        error("$(var_name) decision variables appear to be missing for $(comp_name) components")
    end
end


"""
    @smart_constraint model::JuMP.Model vars::Vector expr::JuMP.Expression

Detection of whether a constraint should be NL or not"
"""
macro smart_constraint(model, vars, expr)
    esc(quote
        if _has_nl_expression($vars)
            JuMP.@NLconstraint($model, $expr)
        else
            JuMP.@constraint($model, $expr)
        end
    end)
end


"""
    set_lower_bound(x::JuMP.VariableRef, bound)

Local wrapper method for JuMP.set_lower_bound, which skips NaN and infinite (-Inf only)
"""
function set_lower_bound(x::JuMP.VariableRef, bound; loose_bounds::Bool=false, pm=missing, category::Symbol=:default)
    if !(isnan(bound) || bound==-Inf)
        JuMP.set_lower_bound(x, bound)
    elseif loose_bounds
        lbs = pm.ext[:loose_bounds]
        JuMP.set_lower_bound(x, -lbs.bound_values[category])
        push!(lbs.loose_lb_vars, x)
    end
end


"""
    set_upper_bound(x::JuMP.VariableRef, bound)

Local wrapper method for JuMP.set_upper_bound, which skips NaN and infinite (+Inf only)
"""
function set_upper_bound(x::JuMP.VariableRef, bound; loose_bounds::Bool=false, pm=missing, category::Symbol=:default)
    if !(isnan(bound) || bound==Inf)
        JuMP.set_upper_bound(x, bound)
    elseif loose_bounds
        lbs = pm.ext[:loose_bounds]
        JuMP.set_upper_bound(x, lbs.bound_values[category])
        push!(lbs.loose_ub_vars, x)
    end
end


"""
    comp_start_value(comp::Dict, keys::Vector{String}, conductor::Int, default::Any)

Searches for start value for a variable `key` in order from a list of `keys` of a component `comp`
for conductor `conductor`, and if one does not exist, uses `default`
"""
function comp_start_value(comp::Dict{String,<:Any}, keys::Vector{String}, conductor::Int, default)
    cond_ind = _get_conductor_indicator(comp)
    for key in keys
        if haskey(comp, key) && !isempty(cond_ind)
            return comp[key][findfirst(isequal(conductor), comp[cond_ind])]
        end
    end
    return default
end


"""
    comp_start_value(comp::Dict, key::String, conductor::Int, default::Any)

Searches for start value for a variable `key` of a component `comp` for conductor `conductor`,
and if one does not exist, uses `default`
"""
function comp_start_value(comp::Dict{String,<:Any}, key::String, conductor::Int, default)
    cond_ind = _get_conductor_indicator(comp)
    if haskey(comp, key) && !isempty(cond_ind)
        return comp[key][findfirst(isequal(conductor), comp[cond_ind])]
    else
        return default
    end
end


"""
    comp_start_value(comp::Dict, keys::Vector{String}, default::Any)

Searches for start value for a variable `key` in order from a list of `keys` of a component `comp`,
and if one does not exist, uses `default`. This is the conductor-agnostic version of `comp_start_value`.
"""
function comp_start_value(comp::Dict{String,<:Any}, keys::Vector{String}, default=0.0)
    for key in keys
        if haskey(comp, key)
            return comp[key]
        end
    end
    return default
end


"""
    comp_start_value(comp::Dict, key::String, default::Any)

Searches for start value for a variable `key` of a component `comp`, and if one does not exist,
uses `default`. This is the conductor-agnostic version of `comp_start_value`.
"""
function comp_start_value(comp::Dict{String,<:Any}, key::String, default=0.0)
    return get(comp, key, default)
end
