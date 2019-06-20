export LPLinUBFPowerModel, LPLinUBFForm

"LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current in scalar form"
abstract type LPLinUBFForm <: _PMs.AbstractBFForm end

""
const LPLinUBFPowerModel = _PMs.GenericPowerModel{LPLinUBFForm}

"default Lin3Distflow constructor for scalar form"
LPLinUBFPowerModel(data::Dict{String,Any}; kwargs...) =
    _PMs.GenericPowerModel(data, LPLinUBFForm; kwargs...)


"rolls a 1d array left or right by idx"
function roll(array::Array{T, 1}, idx::Int; right=true) where T <: Number
    out = Array{T}(undef, size(array))
    pos = idx % length(out)

    if right
        out[1+pos:end] = array[1:end-pos]
        out[1:pos] = array[end-(pos-1):end]
    else
        out[1:end-pos] = array[1+pos:end]
        out[end-(pos-1):end] = array[1:pos]
    end

    return out
end


""
function variable_tp_voltage(pm::_PMs.GenericPowerModel{T}; kwargs...) where T <: LPLinUBFForm
    for cnd in _PMs.conductor_ids(pm)
        _PMs.variable_voltage_magnitude_sqr(pm; cnd=cnd, kwargs...)
    end
end

function variable_tp_branch_current(pm::_PMs.GenericPowerModel{T}; kwargs...) where T <: LPLinUBFForm
    # nothing to do, variables not used in linearised branch flow model
end


"""
Defines voltage drop over a branch, linking from and to side voltage magnitude
"""
function constraint_tp_voltage_magnitude_difference(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm) where T <: LPLinUBFForm
    p_fr = [_PMs.var(pm, n, d, :p, f_idx) for d in _PMs.conductor_ids(pm)]
    q_fr = [_PMs.var(pm, n, d, :q, f_idx) for d in _PMs.conductor_ids(pm)]
    w_fr = _PMs.var(pm, n, c, :w, f_bus)
    w_to = _PMs.var(pm, n, c, :w, t_bus)

    np = length(_PMs.conductor_ids(pm))
    rot = roll([_wrap_to_pi(2*pi/np*(1-d)) for d in _PMs.conductor_ids(pm)], c-1)

    #KVL over the line:
    JuMP.@constraint(pm.model, w_to == w_fr - 2*sum((r[c,d]*cos(rot[d])-x[c,d]*sin(rot[d]))*(p_fr[d] - g_sh_fr*(w_fr/tm^2)) +
                                               (r[c,d]*sin(rot[d])+x[c,d]*cos(rot[d]))*(q_fr[d] + b_sh_fr*(w_fr/tm^2)) for d in _PMs.conductor_ids(pm)) )
end

""
function _PMs.constraint_model_current(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, f_idx, g_sh_fr, b_sh_fr, tm) where T <: LPLinUBFForm
end


"""
Defines branch flow model power flow loss equations
"""
function constraint_tp_flow_losses(pm::_PMs.GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to) where T <: LPLinUBFForm
    tm = [1,1,1] #TODO
    for c in _PMs.conductor_ids(pm)
        constraint_tp_flow_losses(pm, n, c, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr[c,c], g_sh_to[c,c], b_sh_fr[c,c], b_sh_to[c,c], tm[c])
    end
end



"""
Defines branch flow model power flow equations
"""
function constraint_tp_model_voltage_magnitude_difference(pm::_PMs.GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm) where T <: LPLinUBFForm
    tm = [1,1,1] #TODO
    for c in _PMs.conductor_ids(pm)
        constraint_tp_voltage_magnitude_difference(pm, n, c, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr[c,c], b_sh_fr[c,c], tm[c])
    end
end



"""
Defines branch flow model power flow equations
"""
function constraint_tp_flow_losses(pm::_PMs.GenericPowerModel{T}, n::Int, c::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to, tm) where T <: LPLinUBFForm
    p_fr = _PMs.var(pm, n, c, :p, f_idx)
    q_fr = _PMs.var(pm, n, c, :q, f_idx)
    p_to = _PMs.var(pm, n, c, :p, t_idx)
    q_to = _PMs.var(pm, n, c, :q, t_idx)
    w_fr = _PMs.var(pm, n, c, :w, f_bus)
    w_to = _PMs.var(pm, n, c, :w, t_bus)

    JuMP.@constraint(pm.model, p_fr + p_to ==  g_sh_fr*(w_fr/tm^2) +  g_sh_to*w_to)
    JuMP.@constraint(pm.model, q_fr + q_to == -b_sh_fr*(w_fr/tm^2) + -b_sh_to*w_to)
end
