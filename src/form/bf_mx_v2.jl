#TODO
#- Fix matrix variable dict keys once merging
#- Add formulation types
#- Add matrix KCL
import LinearAlgebra: diag, diagm


# MATRIX POWER BALANCE
function variable_tp_generation_power_mx(pm::_PMs.GenericPowerModel; nw=pm.cnw)
    gen_ids = collect(_PMs.ids(pm, nw, :gen))
    ncnds = length(_PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(gen_ids), Array{Real,2}}()
    for (id, gen) in _PMs.ref(pm, nw, :gen)
        bus_id = gen["gen_bus"]
        vmax = _PMs.ref(pm, nw, :bus, bus_id)["vmax"].values
        vmin = _PMs.ref(pm, nw, :bus, bus_id)["vmin"].values
        pmax = max(abs.(gen["pmax"].values), abs.(gen["pmin"].values))
        qmax = max(abs.(gen["qmax"].values), abs.(gen["qmin"].values))
        smax = sqrt.(pmax.^2 + qmax.^2)
        cmax = smax./vmin
        bound[id] = vmax*cmax'
    end
    # create matrix variables
    (Pg,Qg) = variable_mx_complex(pm, gen_ids, ncnds, bound; name=("Pg", "Qg"), prefix="$nw")
    # save references
    _PMs.var(pm, nw)[:Pg] = Pg
    _PMs.var(pm, nw)[:Qg] = Qg
    for c in 1:ncnds
        _PMs.var(pm, nw, c)[:pg] =Dict([(id, Pg[id][c,c]) for id in gen_ids])
        _PMs.var(pm, nw, c)[:qg] =Dict([(id, Qg[id][c,c]) for id in gen_ids])
    end
end


# MATRIX POWER BALANCE
function variable_tp_generation_current_mx(pm::_PMs.GenericPowerModel; nw=pm.cnw)
    gen_ids = collect(_PMs.ids(pm, nw, :gen))
    ncnds = length(_PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(gen_ids), Array{Real,2}}()
    for (id, gen) in _PMs.ref(pm, nw, :gen)
        bus_id = gen["gen_bus"]
        vmax = _PMs.ref(pm, nw, :bus, bus_id)["vmax"].values
        vmin = _PMs.ref(pm, nw, :bus, bus_id)["vmin"].values
        pmax = max(abs.(gen["pmax"].values), abs.(gen["pmin"].values))
        qmax = max(abs.(gen["qmax"].values), abs.(gen["qmin"].values))
        smax = sqrt.(pmax.^2 + qmax.^2)
        bound[id] = cmax*cmax'
    end
    # create matrix variables
    (Lre,Lim) = variable_mx_hermitian(pm, gen_ids, ncnds, bound; name="Ld", prefix="$nw")
    # save references
    _PMs.var(pm, nw)[:Lgre] = Lre
    _PMs.var(pm, nw)[:Lgre] = Lim
end


function variable_tp_load_power_mx(pm::_PMs.GenericPowerModel; nw=pm.cnw)
    load_ids = collect(_PMs.ids(pm, nw, :load))
    ncnds = length(_PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(load_ids), Array{Real,2}}()
    for (id, load) in _PMs.ref(pm, nw, :load)
        bus_id = load["load_bus"]
        vmax = _PMs.ref(pm, nw, :bus, bus_id)["vmax"].values
        vmin = _PMs.ref(pm, nw, :bus, bus_id)["vmin"].values
        # this presumes constant power, wye loads!
        @assert(load["model"]=="constant_power")
        #TODO extend to other load models
        pmax = abs.(load["pd"].values)
        qmax = abs.(load["qd"].values)
        smax = sqrt.(pmax.^2 + qmax.^2)
        cmax = smax./vmin
        bound[id] = vmax*cmax'
    end
    # create matrix variables
    (Pd,Qd) = variable_mx_complex(pm, load_ids, ncnds, bound; name=("Pd", "Qd"), prefix="$nw")
    # save references
    _PMs.var(pm, nw)[:Pd] = Pd
    _PMs.var(pm, nw)[:Qd] = Qd
    for c in 1:ncnds
        _PMs.var(pm, nw, c)[:pd] =Dict([(id, Pd[id][c,c]) for id in gen_ids])
        _PMs.var(pm, nw, c)[:qd] =Dict([(id, Qd[id][c,c]) for id in gen_ids])
    end
end


function variable_tp_load_current_mx(pm::_PMs.GenericPowerModel; nw=pm.cnw)
    load_ids = collect(_PMs.ids(pm, nw, :load))
    ncnds = length(_PMs.conductor_ids(pm, nw))
    # calculate bounds
    bound = Dict{eltype(load_ids), Array{Real,2}}()
    for (id, load) in _PMs.ref(pm, nw, :load)
        bus_id = load["load_bus"]
        vmax = _PMs.ref(pm, nw, :bus, bus_id)["vmax"].values
        vmin = _PMs.ref(pm, nw, :bus, bus_id)["vmin"].values
        # this presumes constant power, wye loads!
        @assert(load["model"]=="constant_power")
        #TODO extend to other load models
        pmax = abs.(load["pd"].values)
        qmax = abs.(load["qd"].values)
        smax = sqrt.(pmax.^2 + qmax.^2)
        cmax = smax./vmin
        bound[id] = cmax*cmax'
    end
    # create matrix variables
    (Ldre, Ldim) = variable_mx_hermitian(pm, load_ids, ncnds, bound; name="Ld", prefix="$nw")
    # save references
    _PMs.var(pm, nw)[:Ldre] = Ldre
    _PMs.var(pm, nw)[:Ldim] = Ldim
end


function constraint_tp_generation_current(pm::_PMs.GenericPowerModel{T}, gen_id::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    Pg = _PMs.var(pm, nw, :Pg, gen_id)
    Qg = _PMs.var(pm, nw, :Qg, gen_id)
    bus_id = _PMs.ref(pm, nw, :gen)["gen_bus"]
    W_re = _PMs.var(pm, nw, :W_re, bus_id)
    W_im = _PMs.var(pm, nw, :W_im, bus_id)
    Ldre = _PMs.var(pm, nw, :Lgre, load_id)
    Ldim = _PMs.var(pm, nw, :Lgim, load_id)
    constraint_SWL_psd(Pg, Qg, W_re, W_im, Lgre, Lgim)
end


function constraint_tp_load_current(pm::_PMs.GenericPowerModel{T}, load_id::Int; nw::Int=pm.cnw) where T <: AbstractUBFForm
    Pd = _PMs.var(pm, nw, :Pd, load_id)
    Qd = _PMs.var(pm, nw, :Qd, load_id)
    bus_id = _PMs.ref(pm, nw, :load)["load_bus"]
    W_re = _PMs.var(pm, nw, :W_re, bus_id)
    W_im = _PMs.var(pm, nw, :W_im, bus_id)
    Ldre = _PMs.var(pm, nw, :Ldre, load_id)
    Ldim = _PMs.var(pm, nw, :Ldim, load_id)
    constraint_SWL_psd(pm, Pd, Qd, W_re, W_im, Ldre, Ldim)
end


function constraint_tp_voltage_psd(pm::_PMs.GenericPowerModel; nw=pm.cnw)
    buses_covered = [i for (l,i,j) in _PMs.ref(pm, nw, :arcs)]
    buses_psd = [i for i in _PMs.ids(pm, nw, :bus) if !(i in buses_covered)]
    for bus_id in buses_psd
        W_re = _PMs.var(pm, nw, :W_re, bus_id)
        W_im = _PMs.var(pm, nw, :W_im, bus_id)
        constraint_M_psd(W_re, W_im)
    end
end


function constraint_SWL_psd(pm::_PMs.GenericPowerModel, P, Q, W_re, W_im, L_re, L_im)
    M_re = [W_re P; P' L_re]
    M_im = [W_im Q; -Q' L_im]
    constraint_M_psd(pm, M_re, M_im)
end


function constraint_M_psd(pm, M_re, M_im)
    JuMP.@constraint(pm.model, [M_re -M_im; M_im M_re] in JuMP.PSDCone())
end
