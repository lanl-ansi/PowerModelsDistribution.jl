""
function run_ac_mc_pf_lm(data, solver; kwargs...)
    return run_mc_pf_lm(data, _PMs.ACPPowerModel, solver; kwargs...)
end


""
function run_dc_mc_opf_lm(data, solver; kwargs...)
    return run_mc_pf_lm(data, _PMs.DCPPowerModel, solver; kwargs...)
end


""
function run_mc_pf_lm(data::Dict{String,Any}, model_type, solver; kwargs...)
    return _PMs.run_model(data, model_type, solver, post_mc_pf_lm; multiconductor=true, ref_extensions=[ref_add_arcs_trans!], kwargs...)
end


""
function run_mc_pf_lm(file::String, model_type, solver; kwargs...)
    return run_mc_pf_lm(PowerModelsDistribution.parse_file(file), model_type, solver; kwargs...)
end


""
function post_mc_pf_lm(pm::_PMs.AbstractPowerModel)
    variable_mc_voltage(pm, bounded=false)
    variable_mc_branch_flow(pm, bounded=false)

    for c in _PMs.conductor_ids(pm)
        _PMs.variable_generation(pm, bounded=false, cnd=c)
        variable_mc_load(pm, bounded=false, cnd=c)
        _PMs.variable_dcline_flow(pm, bounded=false, cnd=c)
    end

    variable_mc_trans_flow(pm, bounded=false)

    constraint_mc_model_voltage(pm)

    for (i,bus) in _PMs.ref(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)

        for c in _PMs.conductor_ids(pm)
            @assert bus["bus_type"] == 3
            _PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c)
        end
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in _PMs.ids(pm, :load)
        constraint_mc_load(pm, id)
    end

    for (i,bus) in _PMs.ref(pm, :bus), c in _PMs.conductor_ids(pm)
        constraint_mc_power_balance_load(pm, i, cnd=c)

        # PV Bus Constraints
        if length(_PMs.ref(pm, :bus_gens, i)) > 0 && !(i in _PMs.ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            _PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c)
            for j in _PMs.ref(pm, :bus_gens, i)
                _PMs.constraint_active_gen_setpoint(pm, j, cnd=c)
            end
        end
    end

    for i in _PMs.ids(pm, :branch), c in _PMs.conductor_ids(pm)
        constraint_mc_ohms_yt_from(pm, i, cnd=c)
        constraint_mc_ohms_yt_to(pm, i, cnd=c)
    end

    for (i,dcline) in _PMs.ref(pm, :dcline), c in _PMs.conductor_ids(pm)

        _PMs.constraint_active_dcline_setpoint(pm, i, cnd=c)

        f_bus = _PMs.ref(pm, :bus)[dcline["f_bus"]]
        if f_bus["bus_type"] == 1
            _PMs.constraint_voltage_magnitude_setpoint(pm, f_bus["index"], cnd=c)
        end

        t_bus = _PMs.ref(pm, :bus)[dcline["t_bus"]]
        if t_bus["bus_type"] == 1
            _PMs.constraint_voltage_magnitude_setpoint(pm, t_bus["index"], cnd=c)
        end
    end

    for i in _PMs.ids(pm, :trans)
        constraint_mc_trans(pm, i)
    end
end
