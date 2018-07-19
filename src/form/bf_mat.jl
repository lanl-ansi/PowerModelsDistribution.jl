export
    AbstractConicUBFForm, AbstractNLPUBFForm,
    SDPUBFPowerModel, SDPUBFForm,
    LPUBFForm, LPUBFPowerModel, LPfullUBFForm, LPdiagUBFPowerModel, LPdiagUBFForm,
    SOCUBFForm, SOCConicUBFPowerModel, SOCNLPUBFPowerModel, SOCConicUBFForm, SOCNLPUBFForm

""
abstract type AbstractNLPUBFForm <: PMs.AbstractPowerFormulation end

""
abstract type AbstractConicUBFForm <: PMs.AbstractConicPowerFormulation end

AbstractUBFForm = Union{AbstractNLPUBFForm, AbstractConicUBFForm}

"SDP BFM per Gan and Low 2014, PSCC"
abstract type SDPUBFForm <: AbstractConicUBFForm end


"SOC relaxation of SDPUBFForm per Kim, Kojima, & Yamashita 2003, cast as an NLP"
abstract type SOCNLPUBFForm <: AbstractNLPUBFForm end

"SOC relaxation of SDPUBFForm per Kim, Kojima, & Yamashita 2003, cast as a SOC"
abstract type SOCConicUBFForm <: AbstractConicUBFForm end

SOCUBFForm = Union{SOCNLPUBFForm, SOCConicUBFForm}


"Simplified BFM per Gan and Low 2014, PSCC, using matrix variables for power, voltage and current"
abstract type LPfullUBFForm <: AbstractNLPUBFForm end

"LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current"
abstract type LPdiagUBFForm <: AbstractNLPUBFForm end

LPUBFForm = Union{LPfullUBFForm, LPdiagUBFForm}

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

function variable_tp_branch_current(pm::GenericPowerModel{T}; kwargs...) where T <: AbstractUBFForm
    variable_tp_branch_series_current_prod_hermitian(pm; kwargs...)
end

function variable_tp_voltage(pm::GenericPowerModel{T}; kwargs...) where T <: AbstractUBFForm
    variable_tp_voltage_prod_hermitian(pm; kwargs...)
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
            var(pm, nw, c)[:cm] = @variable(pm.model,
                [l in ids(pm, nw, :branch)], basename="$(nw)_$(c)_cm",
                lowerbound = 0,
                upperbound = (cmax[l])^2,
                start = PMs.getval(ref(pm, nw, :branch, l), "i_start", c) #TODO shouldn't this be squared?
            )
        else
            var(pm, nw, c)[:cm] = @variable(pm.model,
                [l in ids(pm, nw, :branch)], basename="$(nw)_$(c)_cm",
                lowerbound = 0,
                start = PMs.getval(ref(pm, nw, :branch, l), "i_start", c)
            )
        end
        # PowerModels.variable_current_magnitude_sqr(pm, nw=nw, cnd=c)
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
        ccm =  [var(pm, nw, h, :cm,  i) for h in 1:n_diag_el]
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
function variable_tp_branch_flow(pm::GenericPowerModel; n_cond::Int=3, nw::Int=pm.cnw, bounded = true)
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
function constraint_tp_flow_losses(pm::GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to) where T <: AbstractUBFForm
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
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, n::Int, i) where T <: AbstractUBFForm
    nconductors = length(PMs.conductor_ids(pm))

    w_re = var(pm, n, :w_re)[i]
    w_im = var(pm, n, :w_im)[i]

    alpha = exp(-im*ThreePhasePowerModels.wraptopi(2 * pi / nconductors ))
    beta = (alpha*ones(nconductors)).^(0:nconductors-1)
    gamma = beta*beta'

    w_re_ref = real(gamma).*w_re[1,1]
    w_im_ref = imag(gamma).*w_re[1,1]
    @constraint(pm.model, diag(w_re)[2:3]        .== diag(w_re_ref)[2:3]) # first equality is implied
    @constraint(pm.model, mat2utrivec(w_re) .== mat2utrivec(w_re_ref))
    @constraint(pm.model, mat2utrivec(w_im) .== mat2utrivec(w_im_ref))
end

"""
```
sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd[d] for d in bus_loads) - sum(gs[s] for s in bus_shunts)*v^2
sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd[d] for d in bus_loads) + sum(bs[s] for s in bus_shunts)*v^2
```
"""
function constraint_tp_kcl_shunt(pm::GenericPowerModel{T}, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: AbstractUBFForm
    for c in PMs.conductor_ids(pm)
        w = [var(pm, n, c, :w,  i)]

        p = [var(pm, n, c, :p, a) for a in bus_arcs]
        q = [var(pm, n, c, :q, a) for a in bus_arcs]
        pg = [var(pm, n, c, :pg, g) for g in bus_gens]
        qg = [var(pm, n, c, :qg, g) for g in bus_gens]
        p_dc = [var(pm, n, c, :p_dc, d) for d in bus_arcs_dc]
        q_dc = [var(pm, n, c, :q_dc, d) for d in bus_arcs_dc]

        PMs.con(pm, n, c, :kcl_p)[i] = @constraint(pm.model, sum(pac for (i, pac) in enumerate(p)) + sum(pdc for (i, pdc) in enumerate(p_dc)) == sum(pgen for (i, pgen) in enumerate(pg)) - sum(pd[c] for pd in values(bus_pd)) - sum(gs[c] for gs in values(bus_gs))*w)
        PMs.con(pm, n, c, :kcl_q)[i] = @constraint(pm.model, sum(qac for (i, qac) in enumerate(q)) + sum(qdc for (i, qdc) in enumerate(q_dc)) == sum(qgen for (i, qgen) in enumerate(qg)) - sum(qd[c] for qd in values(bus_qd)) + sum(bs[c] for bs in values(bus_bs))*w)
    end
end


"""
Defines voltage drop over a branch, linking from and to side voltage
"""
function constraint_tp_voltage_magnitude_difference(pm::GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, b_sh_fr, tm) where T <: AbstractUBFForm
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

""
function get_solution_tp(pm::GenericPowerModel, sol::Dict{String,Any})
    add_bus_voltage_setpoint(sol, pm)
    PMs.add_generator_power_setpoint(sol, pm)
    add_branch_flow_setpoint(sol, pm)
    PMs.add_dcline_flow_setpoint(sol, pm)

    PMs.add_kcl_duals(sol, pm)
    PMs.add_sm_duals(sol, pm) # Adds the duals of the transmission lines' thermal limits.
end


""
function add_bus_voltage_setpoint(sol, pm::GenericPowerModel)
    PMs.add_setpoint(sol, pm, "bus", "vm", :w; scale = (x,item) -> sqrt(x))
    PMs.add_setpoint(sol, pm, "bus", "w",  :w)
    PMs.add_setpoint(sol, pm, "bus", "wr", :wr)
    PMs.add_setpoint(sol, pm, "bus", "wi", :wi)
end

""
function add_branch_flow_setpoint(sol, pm::GenericPowerModel)
    if haskey(pm.setting, "output") && haskey(pm.setting["output"], "branch_flows") && pm.setting["output"]["branch_flows"] == true
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

        PMs.add_setpoint(sol, pm, "branch", "cm", :cm; scale = (x,item) -> sqrt(x))
        PMs.add_setpoint(sol, pm, "branch", "ccm", :cm)
        PMs.add_setpoint(sol, pm, "branch", "ccmr", :ccmr)
        PMs.add_setpoint(sol, pm, "branch", "ccmi", :ccmi)
    end
end








## Functions below for ompatibility as unbalanced W-space formulations are not yet <: AbstractWForms

""
function constraint_kcl_shunt(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    if !haskey(con(pm, nw, cnd), :kcl_p)
        con(pm, nw, cnd)[:kcl_p] = Dict{Int,ConstraintRef}()
    end
    if !haskey(con(pm, nw, cnd), :kcl_q)
        con(pm, nw, cnd)[:kcl_q] = Dict{Int,ConstraintRef}()
    end

    bus = ref(pm, nw, :bus, i)
    bus_arcs = ref(pm, nw, :bus_arcs, i)
    bus_arcs_dc = ref(pm, nw, :bus_arcs_dc, i)
    bus_gens = ref(pm, nw, :bus_gens, i)
    bus_loads = ref(pm, nw, :bus_loads, i)
    bus_shunts = ref(pm, nw, :bus_shunts, i)

    bus_pd = Dict(k => ref(pm, nw, :load, k, "pd", cnd) for k in bus_loads)
    bus_qd = Dict(k => ref(pm, nw, :load, k, "qd", cnd) for k in bus_loads)

    bus_gs = Dict(k => ref(pm, nw, :shunt, k, "gs", cnd) for k in bus_shunts)
    bus_bs = Dict(k => ref(pm, nw, :shunt, k, "bs", cnd) for k in bus_shunts)

    constraint_kcl_shunt(pm, nw, cnd, i, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs)
end
"""
```
sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd[d] for d in bus_loads) - sum(gs[s] for d in bus_shunts)*w[i]
sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd[d] for d in bus_loads) + sum(bs[s] for d in bus_shunts)*w[i]
```
"""
function constraint_kcl_shunt(pm::GenericPowerModel{T}, n::Int, c::Int, i, bus_arcs, bus_arcs_dc, bus_gens, bus_pd, bus_qd, bus_gs, bus_bs) where T <: AbstractUBFForm
    w    = var(pm, n, c, :w, i)
    pg   = var(pm, n, c, :pg)
    qg   = var(pm, n, c, :qg)
    p    = var(pm, n, c, :p)
    q    = var(pm, n, c, :q)
    p_dc = var(pm, n, c, :p_dc)
    q_dc = var(pm, n, c, :q_dc)

    @constraint(pm.model, sum(p[a] for a in bus_arcs) + sum(p_dc[a_dc] for a_dc in bus_arcs_dc) == sum(pg[g] for g in bus_gens) - sum(pd for pd in values(bus_pd)) - sum(gs for gs in values(bus_gs))*w)
    @constraint(pm.model, sum(q[a] for a in bus_arcs) + sum(q_dc[a_dc] for a_dc in bus_arcs_dc) == sum(qg[g] for g in bus_gens) - sum(qd for qd in values(bus_qd)) + sum(bs for bs in values(bus_bs))*w)
end


""
function constraint_voltage_angle_difference(pm::GenericPowerModel, i::Int; nw::Int=pm.cnw, cnd::Int=pm.ccnd)
    branch = ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    pair = (f_bus, t_bus)
    buspair = ref(pm, nw, :buspairs, pair)

    if buspair["branch"] == i
        constraint_voltage_angle_difference(pm, nw, cnd, f_idx, buspair["angmin"][cnd], buspair["angmax"][cnd])
    end
end


function constraint_voltage_angle_difference(pm::GenericPowerModel{T}, n::Int, c::Int, f_idx, angmin, angmax) where T <: AbstractUBFForm
    i, f_bus, t_bus = f_idx
    t_idx = (i, t_bus, f_bus)

    branch = ref(pm, n, :branch, i)
    tm = branch["tap"][c]
    g, b = PowerModels.calc_branch_y(branch)
    g, b = g[c,c], b[c,c]
    g_fr = branch["g_fr"][c]
    g_to = branch["g_to"][c]
    b_fr = branch["b_fr"][c]
    b_to = branch["b_to"][c]

    tr, ti = PowerModels.calc_branch_t(branch)
    tr, ti = tr[c], ti[c]

    # convert series admittance to impedance
    z = 1/(g + im*b)
    r = real(z)
    x = imag(z)

    # getting the variables
    w_fr = var(pm, n, c, :w, f_bus)
    p_fr = var(pm, n, c, :p, f_idx)
    q_fr = var(pm, n, c, :q, f_idx)

    tzr = r*tr + x*ti
    tzi = r*ti - x*tr

    @constraint(pm.model,
        tan(angmin)*((tr + tzr*g_fr + tzi*b_fr)*(w_fr/tm^2) - tzr*p_fr + tzi*q_fr)
                 <= ((ti + tzi*g_fr - tzr*b_fr)*(w_fr/tm^2) - tzi*p_fr - tzr*q_fr)
        )
    @constraint(pm.model,
        tan(angmax)*((tr + tzr*g_fr + tzi*b_fr)*(w_fr/tm^2) - tzr*p_fr + tzi*q_fr)
                 >= ((ti + tzi*g_fr - tzr*b_fr)*(w_fr/tm^2) - tzi*p_fr - tzr*q_fr)
        )
end
