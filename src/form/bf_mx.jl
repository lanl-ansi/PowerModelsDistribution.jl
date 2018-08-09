export
    SDPUBFPowerModel, SDPUBFForm,
    LPUBFForm, LPfullUBFPowerModel, LPfullUBFForm, LPdiagUBFPowerModel, LPdiagUBFForm,
    SOCUBFForm, SOCConicUBFPowerModel, SOCNLPUBFPowerModel, SOCConicUBFForm, SOCNLPUBFForm

""
abstract type AbstractNLPUBFForm <: PMs.AbstractBFQPForm end

""
abstract type AbstractConicUBFForm <: PMs.AbstractBFConicForm end

AbstractUBFForm = Union{AbstractNLPUBFForm, AbstractConicUBFForm}


"SDP BFM per Gan and Low 2014, PSCC"
abstract type SDPUBFForm <: AbstractConicUBFForm end


"SOC relaxation of SDPUBFForm per Kim, Kojima, & Yamashita 2003, cast as an QCP"
abstract type SOCNLPUBFForm <: AbstractNLPUBFForm end

"SOC relaxation of SDPUBFForm per Kim, Kojima, & Yamashita 2003, cast as a SOC"
abstract type SOCConicUBFForm <: AbstractConicUBFForm end

SOCUBFForm = Union{SOCNLPUBFForm, SOCConicUBFForm}


"Abstract form for linear unbalanced power flow models"
abstract type AbstractLPUBFForm <: AbstractNLPUBFForm end

"Simplified BFM per Gan and Low 2014, PSCC, using matrix variables for power, voltage and current"
abstract type LPfullUBFForm <: AbstractLPUBFForm end

"LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current"
abstract type LPdiagUBFForm <: AbstractLPUBFForm end


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
const LPfullUBFPowerModel = GenericPowerModel{LPfullUBFForm}

"default LP unbalanced DistFlow constructor"
LPfullUBFPowerModel(data::Dict{String,Any}; kwargs...) =
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
        PowerModels.variable_voltage_magnitude_sqr(pm, nw=nw, cnd=c, bounded=bounded)
    end

    wmaxdict = Dict{Int64, Any}()
    for i in ids(pm, nw, :bus)
        wmax = ref(pm, nw, :bus, i, "vmax").values*ref(pm, nw, :bus, i, "vmax").values'
        wmaxltri = ThreePhasePowerModels.mat2ltrivec(wmax)
        wmaxdict[i] = wmaxltri
    end

    for c in 1:n_lower_triangle_el
        if bounded
            var(pm, nw, c)[:wr] = @variable(pm.model,
                [i in ids(pm, nw, :bus)], basename="$(nw)_$(c)_wr",
                lowerbound = -wmaxdict[i][c],
                upperbound =  wmaxdict[i][c],
                start = PMs.getval(ref(pm, nw, :bus, i), "w_start", c, 1.001)
            )
            var(pm, nw, c)[:wi] = @variable(pm.model,
                [i in ids(pm, nw, :bus)], basename="$(nw)_$(c)_wi",
                lowerbound = -wmaxdict[i][c],
                upperbound =  wmaxdict[i][c],
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
    var(pm, nw)[:W_re] = w_re_dict
    var(pm, nw)[:W_im] = w_im_dict
end

function variable_tp_branch_series_current_prod_hermitian(pm::GenericPowerModel{T}; n_cond::Int=3, nw::Int=pm.cnw, bounded = true) where T <: AbstractUBFForm
    branches = ref(pm, nw, :branch)
    buses = ref(pm, nw, :bus)

    n_diag_el = n_cond
    n_lower_triangle_el = Int((n_cond^2 - n_cond)/2)

    cmax = Dict([(key, zeros(n_cond)) for key in keys(branches)])

    for (key, branch) in branches
        bus_fr = buses[branch["f_bus"]]
        bus_to = buses[branch["t_bus"]]

        vmin_fr = bus_fr["vmin"].values
        vmin_to = bus_to["vmin"].values

        vmax_fr = bus_fr["vmax"].values
        vmax_to = bus_to["vmax"].values

        y_fr_mag = abs.(branch["g_fr"].values + im* branch["b_fr"].values)
        y_to_mag = abs.(branch["g_to"].values + im* branch["b_to"].values)

        smax = branch["rate_a"].values
        cmaxfr = smax./vmin_fr + vmax_fr.*y_fr_mag
        cmaxto = smax./vmin_to + vmax_to.*y_to_mag

        cmax[key] = max.(cmaxfr, cmaxto)
    end


    for c in 1:n_diag_el
        if bounded
            var(pm, nw, c)[:cm] = @variable(pm.model,
                [l in ids(pm, nw, :branch)], basename="$(nw)_$(c)_cm",
                lowerbound = 0,
                upperbound = (cmax[l][c])^2,
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
                lowerbound = -(cmax[l][c])^2,
                upperbound = (cmax[l][c])^2,
                start = PMs.getval(ref(pm, nw, :branch, l), "i_start", c) #TODO shouldn't this be squared?
            )
            var(pm, nw, c)[:ccmi] = @variable(pm.model,
                [l in ids(pm, nw, :branch)], basename="$(nw)_$(c)_ccmi",
                lowerbound = -(cmax[l][c])^2,
                upperbound = (cmax[l][c])^2,
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
    var(pm, nw)[:CC_re] = ccm_re_dict
    var(pm, nw)[:CC_im] = ccm_im_dict
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
    var(pm, nw)[:P_mx] = p_mat_dict
    var(pm, nw)[:Q_mx] = q_mat_dict
end


"variable: `p_lt[l,i,j]` for `(l,i,j)` in `arcs`"
function variable_lower_triangle_active_branch_flow(pm::GenericPowerModel; nw::Int=pm.cnw, cnd::Int=pm.ccnd, bounded = true)
    smaxdict = Dict{Tuple{Int64, Int64, Int64}, Any}()

    for (l,i,j) in ref(pm, nw, :arcs)
        cmax = ref(pm, nw, :branch, l, "rate_a").values./ref(pm, nw, :bus, i, "vmin").values
        smax = ref(pm, nw, :bus, i, "vmax").values.*cmax'
        smaxltri = ThreePhasePowerModels.mat2ltrivec(smax)
        smaxdict[(l,i,j)] = smaxltri
    end

    if bounded
        var(pm, nw, cnd)[:p_lt] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_p_lt",
            lowerbound = -smaxdict[(l,i,j)][cnd],
            upperbound =  smaxdict[(l,i,j)][cnd],
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
    smaxdict = Dict{Tuple{Int64, Int64, Int64}, Any}()

    for (l,i,j) in ref(pm, nw, :arcs)
        cmax = ref(pm, nw, :branch, l, "rate_a").values./ref(pm, nw, :bus, i, "vmin").values
        smax = ref(pm, nw, :bus, i, "vmax").values.*cmax'
        smaxltri = ThreePhasePowerModels.mat2ltrivec(smax)
        smaxdict[(l,i,j)] = smaxltri
    end

    if bounded
        var(pm, nw, cnd)[:q_lt] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_q_lt",
            lowerbound = -smaxdict[(l,i,j)][cnd],
            upperbound =  smaxdict[(l,i,j)][cnd],
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
    smaxdict = Dict{Tuple{Int64, Int64, Int64}, Any}()

    for (l,i,j) in ref(pm, nw, :arcs)
        cmax = ref(pm, nw, :branch, l, "rate_a").values./ref(pm, nw, :bus, i, "vmin").values
        smax = ref(pm, nw, :bus, i, "vmax").values.*cmax'
        smaxutri = ThreePhasePowerModels.mat2utrivec(smax)
        smaxdict[(l,i,j)] = smaxutri
    end

    if bounded
        var(pm, nw, cnd)[:p_ut] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_p_ut",
            lowerbound = -smaxdict[(l,i,j)][cnd],
            upperbound =  smaxdict[(l,i,j)][cnd],
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
    smaxdict = Dict{Tuple{Int64, Int64, Int64}, Any}()

    for (l,i,j) in ref(pm, nw, :arcs)
        cmax = ref(pm, nw, :branch, l, "rate_a").values./ref(pm, nw, :bus, i, "vmin").values
        smax = ref(pm, nw, :bus, i, "vmax").values.*cmax'
        smaxutri = ThreePhasePowerModels.mat2utrivec(smax)
        smaxdict[(l,i,j)] = smaxutri
    end
    if bounded
        var(pm, nw, cnd)[:q_ut] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_q_ut",
            lowerbound = -smaxdict[(l,i,j)][cnd],
            upperbound =  smaxdict[(l,i,j)][cnd],
            start = PMs.getval(ref(pm, nw, :branch, l), "q_start", cnd)
        )
    else
        var(pm, nw, cnd)[:q_ut] = @variable(pm.model,
            [(l,i,j) in ref(pm, nw, :arcs)], basename="$(nw)_$(cnd)_q_ut",
            start = PMs.getval(ref(pm, nw, :branch, l), "q_start", cnd)
        )
    end
end

""
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, i::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    constraint_tp_theta_ref(pm, nw, i)
end


"""
Defines branch flow model power flow equations
"""
function constraint_tp_flow_losses(pm::GenericPowerModel{T}, n::Int, i, f_bus, t_bus, f_idx, t_idx, r, x, g_sh_fr, g_sh_to, b_sh_fr, b_sh_to) where T <: AbstractUBFForm
    p_to = var(pm, n, :P_mx)[t_idx]
    q_to = var(pm, n, :Q_mx)[t_idx]

    p_fr = var(pm, n, :P_mx)[f_idx]
    q_fr = var(pm, n, :Q_mx)[f_idx]

    w_to_re = var(pm, n, :W_re)[t_bus]
    w_fr_re = var(pm, n, :W_re)[f_bus]

    w_to_im = var(pm, n, :W_im)[t_bus]
    w_fr_im = var(pm, n, :W_im)[f_bus]

    ccm_re =  var(pm, n, :CC_re)[i]
    ccm_im =  var(pm, n, :CC_im)[i]

    @constraint(pm.model, p_fr + p_to .==  w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)' + r*ccm_re - x*ccm_im +  w_to_re*(g_sh_to)'  + w_to_im*(b_sh_to)')
    @constraint(pm.model, q_fr + q_to .==  w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)' + x*ccm_re + r*ccm_im +  w_to_im*(g_sh_to)'  - w_to_re*(b_sh_to)')
end

""
function constraint_tp_theta_ref(pm::GenericPowerModel{T}, n::Int, i) where T <: AbstractUBFForm
    nconductors = length(PMs.conductor_ids(pm))

    w_re = var(pm, n, :W_re)[i]
    w_im = var(pm, n, :W_im)[i]

    alpha = exp(-im*ThreePhasePowerModels.wraptopi(2 * pi / nconductors ))
    beta = (alpha*ones(nconductors)).^(0:nconductors-1)
    gamma = beta*beta'

    w_re_ref = real(gamma).*w_re[1,1]
    w_im_ref = imag(gamma).*w_re[1,1]
    @constraint(pm.model, diag(w_re)[2:nconductors]        .== diag(w_re_ref)[2:nconductors]) # first equality is implied
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
    w_fr_re = var(pm, n, :W_re)[f_bus]
    w_fr_im = var(pm, n, :W_im)[f_bus]

    w_to_re = var(pm, n, :W_re)[t_bus]
    w_to_im = var(pm, n, :W_im)[t_bus]

    p_fr = var(pm, n, :P_mx)[f_idx]
    q_fr = var(pm, n, :Q_mx)[f_idx]

    p_s_fr = p_fr - (w_fr_re*(g_sh_fr)' + w_fr_im*(b_sh_fr)')
    q_s_fr = q_fr - (w_fr_im*(g_sh_fr)' - w_fr_re*(b_sh_fr)')

    ccm_re =  var(pm, n, :CC_re)[i]
    ccm_im =  var(pm, n, :CC_im)[i]

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

    add_original_variables(sol, pm)
end


""
function add_bus_voltage_setpoint(sol, pm::GenericPowerModel)
    # PMs.add_setpoint(sol, pm, "bus", "vm_p", :w; scale = (x,item) -> sqrt(x))
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

        # PMs.add_setpoint(sol, pm, "branch", "cm_p", :cm; scale = (x,item) -> sqrt(x))
        # PMs.add_setpoint(sol, pm, "branch", "ccm", :cm)
        # PMs.add_setpoint(sol, pm, "branch", "ccmr", :ccmr)
        # PMs.add_setpoint(sol, pm, "branch", "ccmi", :ccmi)
    end
end


""
function add_original_variables(sol, pm::GenericPowerModel)
    if haskey(pm.setting, "output") && haskey(pm.setting["output"], "original_variables") && pm.setting["output"]["original_variables"] == true
        if !haskey(pm.setting, "output") || !haskey(pm.setting["output"], "branch_flows") || pm.setting["output"]["branch_flows"] == false
            error(LOGGER, "deriving the original variables requires setting: branch_flows = true")
        end

        for (nw, network) in pm.ref[:nw]
            #find rank-1 starting point
            ref_buses = find_ref_buses(pm, nw)
            #TODO develop code to start with any rank-1 W variable
            buses   = ref(pm, nw, :bus)
            arcs    = ref(pm, nw, :arcs)
            branches    = ref(pm, nw, :branch)
            #define sets to explore
            all_bus_ids             = Set([b for (b, bus)    in ref(pm, nw, :bus)])
            all_arc_from_ids        = Set([(l,i,j) for (l,i,j) in ref(pm, nw, :arcs_from)])
            all_arc_to_ids          = Set([(l,i,j) for (l,i,j) in ref(pm, nw, :arcs_to)])
            all_branch_ids          = Set([l for (l,i,j) in ref(pm, nw, :arcs_from)])
            visited_arc_from_ids    = Set()
            visited_arc_to_ids      = Set()
            visited_bus_ids         = Set()
            visited_branch_ids      = Set()

            for b in ref_buses
                # sol["bus"]["$b"]["va"] = wraptopi(ref(pm, nw, :bus, b)["va"].values)
                sol["bus"]["$b"]["va"] = [0, -2*pi/3, 2*pi/3] #TODO support arbitrary angles at the reference bus

                sol["bus"]["$b"]["vm"] = ref(pm, nw, :bus, b)["vm"].values
                # sol["bus"]["$b"]["v"] =  ref(pm, nw, :bus, b)["vm"].values.* exp.(im*sol["bus"]["$b"]["va"])
                push!(visited_bus_ids, b)
            end

            tt = 0
            while visited_branch_ids != all_branch_ids && visited_bus_ids != all_bus_ids
                tt = tt+1
                if tt >100
                    @show "break while"
                    break
                end

                remaining_arc_from_ids = setdiff(all_arc_from_ids, visited_arc_from_ids)
                remaining_arc_to_ids = setdiff(all_arc_to_ids, visited_arc_to_ids)

                candidate_arcs_from = [(l,i,j) for (l,i,j) in remaining_arc_from_ids if i in visited_bus_ids && !(j in visited_bus_ids)]
                candidate_arcs_to   = [(l,i,j) for (l,i,j) in remaining_arc_to_ids   if i in visited_bus_ids && !(j in visited_bus_ids)]


                if !isempty(candidate_arcs_from)
                    (l,i,j) = arc = candidate_arcs_from[1]
                    g_fr = diagm(branches[l]["g_fr"].values)
                    b_fr = diagm(branches[l]["b_fr"].values)
                    y_fr = g_fr + im* b_fr
                    g_to = diagm(branches[l]["g_to"].values)
                    b_to = diagm(branches[l]["b_to"].values)
                    y_to = g_to + im* b_to
                    r = branches[l]["br_r"].values
                    x = branches[l]["br_x"].values
                    z = (r + im*x)
                    Ui = sol["bus"]["$i"]["vm"].*exp.(im*sol["bus"]["$i"]["va"])

                    Pij = make_full_matrix_variable(sol["branch"]["$l"]["pf"].values, sol["branch"]["$l"]["pf_lt"].values, sol["branch"]["$l"]["pf_ut"].values)
                    Qij = make_full_matrix_variable(sol["branch"]["$l"]["qf"].values, sol["branch"]["$l"]["qf_lt"].values, sol["branch"]["$l"]["qf_ut"].values)
                    Sij = Pij + im*Qij

                    Ssij = Sij - Ui*Ui'*y_fr'
                    Isij = (1/trace(Ui*Ui'))*(Ssij')*Ui
                    Uj = Ui - z*Isij
                    Iij = Isij + y_fr*Ui

                    Isji = -Isij
                    Iji = Isji + y_to*Uj


                    sol["bus"]["$j"]["vm"] = abs.(Uj)
                    sol["bus"]["$j"]["va"] = wraptopi(angle.(Uj))

                    sol["branch"]["$l"]["cfm"] = abs.(Iij)
                    sol["branch"]["$l"]["cfa"] = wraptopi(angle.(Iij))
                    sol["branch"]["$l"]["ctm"] = abs.(Iji)
                    sol["branch"]["$l"]["cta"] = wraptopi(angle.(Iji))
#
                    push!(visited_arc_from_ids, arc)
                    push!(visited_branch_ids, l)
                    push!(visited_bus_ids, j)
                    
                elseif !isempty(candidate_arcs_to)
                    (l,i,j) = arc = candidate_arcs_to[1]
                    g_fr = diagm(branches[l]["g_to"].values)
                    b_fr = diagm(branches[l]["b_to"].values)
                    y_fr = g_fr + im* b_fr
                    g_to = diagm(branches[l]["g_fr"].values)
                    b_to = diagm(branches[l]["b_fr"].values)
                    y_to = g_to + im* b_to
                    r = branches[l]["br_r"].values
                    x = branches[l]["br_x"].values
                    z = (r + im*x)
                    @show l, i, j
                    Ui = sol["bus"]["$i"]["vm"].*exp.(im*sol["bus"]["$i"]["va"])

                    Pij = make_full_matrix_variable(sol["branch"]["$l"]["pt"].values, sol["branch"]["$l"]["pt_lt"].values, sol["branch"]["$l"]["pt_ut"].values)
                    Qij = make_full_matrix_variable(sol["branch"]["$l"]["qt"].values, sol["branch"]["$l"]["qt_lt"].values, sol["branch"]["$l"]["qt_ut"].values)
                    Sij = Pij + im*Qij

                    Ssij = Sij - Ui*Ui'*y_fr'
                    Isij = (1/trace(Ui*Ui'))*(Ssij')*Ui
                    Uj = Ui - z*Isij
                    Iij = Isij + y_fr*Ui

                    Isji = -Isij
                    Iji = Isji + y_to*Uj


                    sol["bus"]["$j"]["vm"] = abs.(Uj)
                    sol["bus"]["$j"]["va"] = wraptopi(angle.(Uj))

                    sol["branch"]["$l"]["ctm"] = abs.(Iij)
                    sol["branch"]["$l"]["cta"] = wraptopi(angle.(Iij))
                    sol["branch"]["$l"]["cfm"] = abs.(Iji)
                    sol["branch"]["$l"]["cfa"] = wraptopi(angle.(Iji))
#
                    push!(visited_arc_to_ids, arc)
                    push!(visited_branch_ids, l)
                    push!(visited_bus_ids, j)

                else
                    candidate_arcs = [(l,i,j) for (l,i,j) in remaining_arc_from_ids if i in visited_bus_ids && j in visited_bus_ids]
                    (l,i,j) = arc = candidate_arcs[1]
                    Sij = sol["branch"]["$l"]["pf"] + im* sol["branch"]["$l"]["qf"]
                    Sji = sol["branch"]["$l"]["pt"] + im* sol["branch"]["$l"]["qt"]
                    Ui = sol["bus"]["$i"]["vm"].*exp.(im*sol["bus"]["$i"]["va"])
                    Uj = sol["bus"]["$j"]["vm"].*exp.(im*sol["bus"]["$j"]["va"])

                    Iij = Sij./Ui
                    Iji = Sji./Uj
                    sol["branch"]["$l"]["cfm"] = abs.(Iij)
                    sol["branch"]["$l"]["cfa"] = wraptopi(angle.(Iij))
                    sol["branch"]["$l"]["ctm"] = abs.(Iji)
                    sol["branch"]["$l"]["cta"] = wraptopi(angle.(Iji))
                    push!(visited_arc_from_ids, arc)
                    push!(visited_arc_to_ids, (l,j,i))
                    push!(visited_branch_ids, l)
                end
            end
        end
    end
end
