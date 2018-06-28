"rolls a 1d array left or right by idx"
function roll(array::Array{T, 1}, idx::Int; right=true) where T <: Number
    out = Array{T}(size(array))
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
function variable_tp_voltage(pm::GenericPowerModel{T}; kwargs...) where T <: PMs.AbstractBFForm
    PMs.variable_voltage_magnitude_sqr(pm; kwargs...)
end


"""
Defines voltage drop over a branch, linking from and to side voltage magnitude
"""
function constraint_tp_voltage_magnitude_difference(pm::GenericPowerModel{T}, n::Int, h::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm) where T <: PMs.AbstractBFForm
    p_fr = [var(pm, n, j, :p, f_idx) for j in PMs.phase_ids(pm)]
    q_fr = [var(pm, n, j, :q, f_idx) for j in PMs.phase_ids(pm)]
    w_fr = var(pm, n, h, :w, f_bus)
    w_to = var(pm, n, h, :w, t_bus)

    np = length(PMs.phase_ids(pm))
    rot = roll([2*pi/np*(j-1) for j in PMs.phase_ids(pm)], h-1)

    #KVL over the line:
    @constraint(pm.model, w_to == w_fr - 2*sum((r[h,j]*cos(rot[j])-x[h,j]*sin(rot[j]))*p_fr[j] +
                                               (r[h,j]*sin(rot[j])+x[h,j]*cos(rot[j]))*q_fr[j] for j in PMs.phase_ids(pm)) )
end
