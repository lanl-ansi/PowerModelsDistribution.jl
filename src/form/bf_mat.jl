export
    AbstractConicUBFForm, AbstractNLPUBFForm,
    SDPUBFPowerModel, SDPUBFForm,
    LPUBFForm, LPUBFPowerModel, LPfullUBFForm, LPdiagUBFPowerModel, LPdiagUBFForm,
    SOCUBFForm, SOCConicUBFPowerModel, SOCNLPUBFPowerModel, SOCConicUBFForm, SOCNLPUBFForm

""
abstract type AbstractNLPUBFForm <: PMs.AbstractPowerFormulation end

""
abstract type AbstractConicUBFForm <: PMs.AbstractConicPowerFormulation end


""
abstract type SOCNLPUBFForm <: AbstractNLPUBFForm end

""
abstract type SOCConicUBFForm <: AbstractConicUBFForm end

SOCUBFForm = Union{SOCNLPUBFForm, SOCConicUBFForm}


""
abstract type LPfullUBFForm <: AbstractNLPUBFForm end

""
abstract type LPdiagUBFForm <: AbstractNLPUBFForm end

LPUBFForm = Union{LPfullUBFForm, LPdiagUBFForm}

""
abstract type SDPUBFForm <: AbstractConicUBFForm end

AbstractUBFForm = Union{AbstractNLPUBFForm, AbstractConicUBFForm}


""
const SDPUBFPowerModel = GenericPowerModel{SDPUBFForm}

"default SDP unbalanced DistFlow constructor"
SDPUBFPowerModel(data::Dict{String,Any}; kwargs...) =
    GenericPowerModel(data, SDPUBFForm; kwargs...)

""
const SOCNLPUBFPowerModel = GenericPowerModel{SOCNLPUBFForm}

"default SOC unbalanced DistFlow constructor"
SOCNLPUBFPowerModel(data::Dict{String,Any}; kwargs...) =
GenericPowerModel(data, SOCNLPUBFForm; kwargs...)

""
const SOCConicUBFPowerModel = GenericPowerModel{SOCConicUBFForm}

"default SOC unbalanced DistFlow constructor"
SOCConicUBFPowerModel(data::Dict{String,Any}; kwargs...) =
GenericPowerModel(data, SOCConicUBFForm; kwargs...)

""
const LPUBFPowerModel = GenericPowerModel{LPfullUBFForm}

"default LP unbalanced DistFlow constructor"
LPUBFPowerModel(data::Dict{String,Any}; kwargs...) =
    GenericPowerModel(data, LPfullUBFForm; kwargs...)

""
const LPdiagUBFPowerModel = GenericPowerModel{LPdiagUBFForm}

"default LP unbalanced DistFlow constructor"
LPdiagUBFPowerModel(data::Dict{String,Any}; kwargs...) =
    GenericPowerModel(data, LPdiagUBFForm; kwargs...)



function variable_branch_current(pm::GenericPowerModel{T}; kwargs...) where T <: AbstractUBFForm
    variable_tp_branch_series_current_prod_hermitian(pm; kwargs...)
end

function variable_tp_voltage(pm::GenericPowerModel{T}; kwargs...) where T <: AbstractUBFForm
    variable_tp_voltage_prod_hermitian(pm; kwargs...)
end

"do nothing, no way to represent this in these variables"
function constraint_tp_voltage_angle_ref(pm::GenericPowerModel{T}, n::Int, h::Int, i, varef) where T <: PMs.AbstractWForms
end


""
function variable_tp_voltage_prod_hermitian(pm::GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    n_diag_el = n_cond
    n_lower_triangle_el = Int((n_cond^2 - n_cond)/2)
    for c in 1:n_diag_el
        PowerModels.variable_voltage_magnitude_sqr(pm, nw=nw, cnd=c)
    end
    for c in 1:n_lower_triangle_el
        if bounded
            var(pm, nw, c)[:wr] = @variable(pm.model,
                [i in ids(pm, nw, :bus)], basename="$(nw)_$(c)_wr",
                lowerbound = -ref(pm, nw, :bus, i, "vmax", c)^2,
                upperbound =  ref(pm, nw, :bus, i, "vmax", c)^2,
                start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
            )
            var(pm, nw, c)[:wi] = @variable(pm.model,
                [i in ids(pm, nw, :bus)], basename="$(nw)_$(c)_wi",
                lowerbound = -ref(pm, nw, :bus, i, "vmax", c)^2,
                upperbound =  ref(pm, nw, :bus, i, "vmax", c)^2,
                start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
            )
        else
            var(pm, nw, c)[:wr] = @variable(pm.model,
                [i in ids(pm, nw, :bus)], basename="$(nw)_$(c)_wr",
                start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
            )
            var(pm, nw, c)[:wi] = @variable(pm.model,
                [i in ids(pm, nw, :bus)], basename="$(nw)_$(c)_wi",
                start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
            )
        end
    end

    #Store dictionary with matrix variables by bus
    w_re_dict = Dict{Int64, Any}()
    w_im_dict = Dict{Int64, Any}()
    for i in ids(pm, nw, :bus)
        w =  [var(pm, nw, h, :w,  i) for h in 1:n_diag_el]
        wr = [var(pm, nw, h, :wr, i) for h in 1:n_lower_triangle_el]
        wi = [var(pm, nw, h, :wi, i) for h in 1:n_lower_triangle_el]

        (w_re, w_im) = make_hermitian_matrix_variable(w, wr, wi)
        w_re_dict[i] = w_re
        w_im_dict[i] = w_im
    end
    var(pm, nw)[:w_re] = w_re_dict
    var(pm, nw)[:w_im] = w_im_dict
end

function variable_tp_branch_series_current_prod_hermitian(pm::GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    branches = ref(pm, nw, :branch)
    buses = ref(pm, nw, :bus)
    assert(n_cond<=5)

    n_diag_el = n_cond
    n_lower_triangle_el = Int((n_cond^2 - n_cond)/2)

    cmax = Dict([(key, 0.0) for key in keys(branches)])

    for (key, branch) in branches
        bus_fr = buses[branch["f_bus"]]
        bus_to = buses[branch["t_bus"]]

        vmin_fr = minimum(bus_fr["vmin"])
        vmin_to = minimum(bus_fr["vmin"])

        smax = maximum(branch["rate_a"])

        M = 2
        cmax[key] = M*smax/min(vmin_fr, vmin_to)
    end


    for c in 1:n_diag_el
        if bounded
            var(pm, nw, c)[:ccm] = @variable(pm.model,
                [l in ids(pm, nw, :branch)], basename="$(nw)_$(c)_ccm",
                lowerbound = 0,
                upperbound = (cmax[l])^2,
                start = PMs.getval(ref(pm, nw, :branch, l), "i_start", c) #TODO shouldn't this be squared?
            )
        else
            var(pm, nw, c)[:ccm] = @variable(pm.model,
                [l in ids(pm, nw, :branch)], basename="$(nw)_$(c)_ccm",
                lowerbound = 0,
                start = PMs.getval(ref(pm, nw, :branch, l), "i_start", c)
            )
        end
    end

    for c in 1:n_lower_triangle_el
        if bounded
            var(pm, nw, c)[:ccmr] = @variable(pm.model,
                [l in ids(pm, nw, :branch)], basename="$(nw)_$(c)_ccmr",
                lowerbound = -(cmax[l])^2,
                upperbound = (cmax[l])^2,
                start = PMs.getval(ref(pm, nw, :branch, l), "i_start", c) #TODO shouldn't this be squared?
            )
            var(pm, nw, c)[:ccmi] = @variable(pm.model,
                [l in ids(pm, nw, :branch)], basename="$(nw)_$(c)_ccmi",
                lowerbound = -(cmax[l])^2,
                upperbound = (cmax[l])^2,
                start = PMs.getval(ref(pm, nw, :branch, l), "i_start", c)
            )
        else
            var(pm, nw, c)[:ccmr] = @variable(pm.model,
                [l in ids(pm, nw, :branch)], basename="$(nw)_$(c)_ccmr",
                # lowerbound = 0,
                start = PMs.getval(ref(pm, nw, :branch, l), "i_start", c)
            )
            var(pm, nw, c)[:ccmi] = @variable(pm.model,
                [l in ids(pm, nw, :branch)], basename="$(nw)_$(c)_ccmi",
                # lowerbound = 0,
                start = PMs.getval(ref(pm, nw, :branch, l), "i_start", c)
            )
        end
    end

    #Store dictionary with matrix variables by branch
    ccm_re_dict = Dict{Int64, Any}()
    ccm_im_dict = Dict{Int64, Any}()
    for i in ids(pm, nw, :branch)
        ccm =  [var(pm, nw, h, :ccm,  i) for h in 1:n_diag_el]
        ccmr = [var(pm, nw, h, :ccmr, i) for h in 1:n_lower_triangle_el]
        ccmi = [var(pm, nw, h, :ccmi, i) for h in 1:n_lower_triangle_el]

        (ccm_re, ccm_im) = make_hermitian_matrix_variable(ccm, ccmr, ccmi)
        ccm_re_dict[i] = ccm_re
        ccm_im_dict[i] = ccm_im
    end
    var(pm, nw)[:ccm_re] = ccm_re_dict
    var(pm, nw)[:ccm_im] = ccm_im_dict
end

""
function variable_branch_flow(pm::GenericPowerModel; n_cond::Int=3, nw::Int=pm.cnw, bounded = true)
    n_diag_el = n_cond
    n_lower_triangle_el = Int((n_cond^2 - n_cond)/2)
    assert(n_cond<=5)

    for i in 1:n_diag_el
        PMs.variable_active_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
        PMs.variable_reactive_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
    end

    for i in 1:n_lower_triangle_el
        variable_lower_triangle_active_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
        variable_lower_triangle_reactive_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
        variable_upper_triangle_active_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
        variable_upper_triangle_reactive_branch_flow(pm, nw=nw, cnd=i, bounded=bounded)
    end
    #Store dictionary with matrix variables by arc
    p_mat_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()
    q_mat_dict = Dict{Tuple{Int64,Int64,Int64}, Any}()

    for i in ref(pm, nw, :arcs)
        p_d =  [var(pm, nw, c, :p,    i) for c in 1:n_diag_el]
        p_ut = [var(pm, nw, c, :p_ut, i) for c in 1:n_lower_triangle_el]
        p_lt = [var(pm, nw, c, :p_lt, i) for c in 1:n_lower_triangle_el]

        q_d =  [var(pm, nw, c, :q,    i) for c in 1:n_diag_el]
        q_ut = [var(pm, nw, c, :q_ut, i) for c in 1:n_lower_triangle_el]
        q_lt = [var(pm, nw, c, :q_lt, i) for c in 1:n_lower_triangle_el]

        p_mat = make_full_matrix_variable(p_d, p_lt, p_ut)
        q_mat = make_full_matrix_variable(q_d, q_lt, q_ut)

        p_mat_dict[i] = p_mat
        q_mat_dict[i] = q_mat
    end
    var(pm, nw)[:p_mat] = p_mat_dict
    var(pm, nw)[:q_mat] = q_mat_dict
end


"variable: `p_lt[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_lower_triangle_active_branch_flow(pm::GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    if bounded
        var(pm, nw, cnd)[:p_lt] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_p_lt",
            lowerbound = -ref(pm, nw, :branch, l, "rate_a", cnd),
            upperbound =  ref(pm, nw, :branch, l, "rate_a", cnd),
            start = PMs.getval(ref(pm, nw, :branch, l), "p_start", cnd)
        )
    else
        var(pm, nw, cnd)[:p_lt] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_p_lt",
            start = PMs.getval(ref(pm, nw, :branch, l), "p_start", cnd)
        )
    end
end

"variable: `q_lt[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_lower_triangle_reactive_branch_flow(pm::GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    if bounded
        var(pm, nw, cnd)[:q_lt] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_q_lt",
            lowerbound = -ref(pm, nw, :branch, l, "rate_a", cnd),
            upperbound =  ref(pm, nw, :branch, l, "rate_a", cnd),
            start = PMs.getval(ref(pm, nw, :branch, l), "q_start", cnd)
        )
    else
        var(pm, nw, cnd)[:q_lt] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_q_lt",
            start = PMs.getval(ref(pm, nw, :branch, l), "q_start", cnd)
        )
    end
end

"variable: `p_ut[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_upper_triangle_active_branch_flow(pm::GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    if bounded
        var(pm, nw, cnd)[:p_ut] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_p_ut",
            lowerbound = -ref(pm, nw, :branch, l, "rate_a", cnd),
            upperbound =  ref(pm, nw, :branch, l, "rate_a", cnd),
            start = PMs.getval(ref(pm, nw, :branch, l), "p_start", cnd)
        )
    else
        var(pm, nw, cnd)[:p_ut] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_p_ut",
            start = PMs.getval(ref(pm, nw, :branch, l), "p_start", cnd)
        )
    end
end

"variable: `q_ut[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_upper_triangle_reactive_branch_flow(pm::GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    if bounded
        var(pm, nw, cnd)[:q_ut] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_q_ut",
            lowerbound = -ref(pm, nw, :branch, l, "rate_a", cnd),
            upperbound =  ref(pm, nw, :branch, l, "rate_a", cnd),
            start = PMs.getval(ref(pm, nw, :branch, l), "q_start", cnd)
        )
    else
        var(pm, nw, cnd)[:q_ut] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_q_ut",
            start = PMs.getval(ref(pm, nw, :branch, l), "q_start", cnd)
        )
    end
end


"""
Defines branch flow model power flow equations
"""
function constraint_tp_flow_losses_mat(pm::GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to) where T <: AbstractUBFForm
    p_to = var(pm, n, :p_mat)[t_idx]
    q_to = var(pm, n, :q_mat)[t_idx]

    p_fr = var(pm, n, :p_mat)[f_idx]
    q_fr = var(pm, n, :q_mat)[f_idx]

    w_to_re = var(pm, n, :w_re)[t_bus]
    w_fr_re = var(pm, n, :w_re)[f_bus]

    ccm_re =  var(pm, n, :ccm_re)[i]
    ccm_im =  var(pm, n, :ccm_im)[i]

    @constraint(pm.model, p_fr + p_to .==  g_sh_fr*(w_fr_re) + r*ccm_re - x*ccm_im +  g_sh_to*w_to_re)
    @constraint(pm.model, q_fr + q_to .== -b_sh_fr*(w_fr_re) + x*ccm_re + r*ccm_im + -b_sh_to*w_to_re)
end

""
function constraint_tp_voltage_ref(pm::GenericPowerModel{T}, n::Int, i, vmref, varef) where T <: AbstractUBFForm
    vref = vmref.*exp.(im*varef)

    w_re = var(pm, n, :w_re)[i]
    w_im = var(pm, n, :w_im)[i]

    wref_re = real(vref*vref')
    wref_im = imag(vref*vref')

    @constraint(pm.model, diag(w_re)        .==        diag(wref_re))
    @constraint(pm.model, mat2utrivec(w_re) .== mat2utrivec(wref_re))
    @constraint(pm.model, mat2utrivec(w_im) .== mat2utrivec(wref_im))
end

"""
```
sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*v^2
sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd[d] for d in bus_loads) + sum(bs[s] for s in bus_shunts)*v^2
```
"""
function constraint_tp_kcl_shunt(pm::GenericPowerModel{T}, n::Int, c::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: AbstractUBFForm
    w = var(pm, n, c, :w, i)
    p = var(pm, n, c, :p)
    q = var(pm, n, c, :q)
    pg = var(pm, n, c, :pg)
    qg = var(pm, n, c, :qg)
    # p_dc = var(pm, n, h, :p_dc)
    # q_dc = var(pm, n, h, :q_dc)

    PMs.con(pm, n, c, :kcl_p)[i] = @constraint(pm.model, sum(p[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*w)
    PMs.con(pm, n, c, :kcl_q)[i] = @constraint(pm.model, sum(q[a] for a in bus_arcs) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*w)
    # PMs.con(pm, n, h, :kcl_p)[i] = @constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*w)
    # PMs.con(pm, n, h, :kcl_q)[i] = @constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*w)
end


"""
Defines voltage drop over a branch, linking from and to side voltage
"""
function constraint_tp_voltage_magnitude_difference_mat(pm::GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm) where T <: AbstractUBFForm
    w_fr_re = var(pm, n, :w_re)[f_bus]
    w_fr_im = var(pm, n, :w_im)[f_bus]

    w_to_re = var(pm, n, :w_re)[t_bus]
    w_to_im = var(pm, n, :w_im)[t_bus]

    p_fr = var(pm, n, :p_mat)[f_idx]
    q_fr = var(pm, n, :q_mat)[f_idx]

    p_s_fr = p_fr - g_sh_fr*w_fr_re
    q_s_fr = q_fr + b_sh_fr*w_fr_re

    ccm_re =  var(pm, n, :ccm_re)[i]
    ccm_im =  var(pm, n, :ccm_im)[i]

        #KVL over the line:
    @constraint(pm.model, diag(w_to_re) .== diag(
                                      w_fr_re   - p_s_fr  *r' - q_s_fr*x'        - r*p_s_fr'    - x*q_s_fr'
                                                + r*ccm_re*r' - x     *ccm_im*r' + x*ccm_re *x' + r*ccm_im *x'))
    @constraint(pm.model, mat2utrivec(w_to_re) .== mat2utrivec(
                                      w_fr_re   - p_s_fr  *r' - q_s_fr*x'        - r*p_s_fr'    - x*q_s_fr'
                                                + r*ccm_re*r' - x     *ccm_im*r' + x*ccm_re *x' + r*ccm_im *x'))
    @constraint(pm.model, mat2utrivec(w_to_im) .== mat2utrivec(
                                      w_fr_im   - q_s_fr  *r' + p_s_fr*x'        - x*p_s_fr'    + r*q_s_fr'
                                                + x*ccm_re*r' + r     *ccm_im*r' - r*ccm_re *x' + x*ccm_im *x'))
end




"""
Defines voltage drop over a branch, linking from and to side voltage
"""
function objective_min_losses(pm::GenericPowerModel{T}) where T <: AbstractUBFForm
    from_idx = Dict()
    to_idx = Dict()
    for (n, nw_ref) in nws(pm)
        from_idx[n] = Dict(arc[1] => arc for arc in nw_ref[:arcs_from])
        to_idx[n] = Dict(arc[1] => arc for arc in nw_ref[:arcs_to])
    end

    return @objective(pm.model, Min,
        sum(
            sum(
                sum(
                    var(pm, n, c, :p, from_idx[n][i]) + var(pm, n, c, :p, to_idx[n][i])
                for c in conductor_ids(pm, n))
            for (i, arc) in nw_ref[:branch])
        for (n, nw_ref) in nws(pm))
    )
end

"""
Defines voltage drop over a branch, linking from and to side voltage
"""
function objective_min_slack(pm::GenericPowerModel)
    slackgen = 1
    return @objective(pm.model, Min,
        sum(
            sum(
                var(pm, n, c, :pg, 1)
            for c in conductor_ids(pm, n))
        for (n, nw_ref) in nws(pm))
    )
end


""
function get_solution_tp(pm::GenericPowerModel, sol::Dict{String,Any})
    add_bus_voltage_setpoint(sol, pm)
    PMs.add_generator_power_setpoint(sol, pm)
    add_branch_flow_setpoint(sol, pm)
    # add_dcline_flow_setpoint(sol, pm)

    PMs.add_kcl_duals(sol, pm)
    PMs.add_sm_duals(sol, pm) # Adds the duals of the transmission lines' thermal limits.
end


""
function add_bus_voltage_setpoint(sol, pm::GenericPowerModel{T}) where T <: AbstractUBFForm
    PMs.add_setpoint(sol, pm, "bus", "vm", :w; scale = (x,item) -> sqrt(x))
    PMs.add_setpoint(sol, pm, "bus", "w",  :w)
    PMs.add_setpoint(sol, pm, "bus", "wr", :wr)
    PMs.add_setpoint(sol, pm, "bus", "wi", :wi)
end

""
function add_branch_flow_setpoint(sol, pm::GenericPowerModel{T}) where T <: AbstractUBFForm
    # check the branch flows were requested
    # if haskey(pm.setting, "output") && haskey(pm.setting["output"], "branch_flows") && pm.setting["output"]["branch_flows"] == true
        PMs.add_setpoint(sol, pm, "branch", "pf", :p; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "qf", :q; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "pf_ut", :p_ut; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "qf_ut", :q_ut; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "pf_lt", :p_lt; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "qf_lt", :q_lt; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])])

        PMs.add_setpoint(sol, pm, "branch", "pt", :p; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "qt", :q; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "pt_ut", :p_ut; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "qt_ut", :q_ut; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "pt_lt", :p_lt; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])])
        PMs.add_setpoint(sol, pm, "branch", "qt_lt", :q_lt; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])])

        PMs.add_setpoint(sol, pm, "branch", "cm", :ccm; scale = (x,item) -> sqrt(x))
        PMs.add_setpoint(sol, pm, "branch", "ccm", :ccm)
        PMs.add_setpoint(sol, pm, "branch", "ccmr", :ccmr)
        PMs.add_setpoint(sol, pm, "branch", "ccmi", :ccmi)

    # end
end
