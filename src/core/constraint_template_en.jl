"""
    function constraint_mc_voltage_absolute(
        pm::RectangularVoltageExplicitNeutralModels,
        id::Int;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true,
        kwargs...
    )

Imposes absolute voltage magnitude bounds for models with explicit neutrals
"""
function constraint_mc_voltage_absolute(pm::RectangularVoltageExplicitNeutralModels, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    bus = ref(pm, nw, :bus, id)

    constraint_mc_voltage_absolute(pm, nw, id, bus["terminals"], bus["grounded"], bus["vmin"], bus["vmax"])
end


"""
    function constraint_mc_voltage_pairwise(
        pm::RectangularVoltageExplicitNeutralModels,
        id::Int;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true,
        kwargs...
    )

Imposes pairwise voltage magnitude bounds, i.e. magnitude bounds on the voltage between to terminals, for models with explicit neutrals
"""
function constraint_mc_voltage_pairwise(pm::RectangularVoltageExplicitNeutralModels, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    bus = ref(pm, nw, :bus, id)

    vm_pair_lb = bus["vm_pair_lb"]
    vm_pair_ub = bus["vm_pair_ub"]

    constraint_mc_voltage_pairwise(pm, nw, id, bus["terminals"], bus["grounded"], vm_pair_lb, vm_pair_ub)
end


"""
    function constraint_mc_voltage_reference(
        pm::ExplicitNeutralModels,
        id::Int;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true,
        kwargs...
    )

Imposes suitable constraints for the voltage at the reference bus
"""
function constraint_mc_voltage_reference(pm::ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    bus = ref(pm, nw, :bus, id)
    terminals = bus["terminals"]
    grounded = bus["grounded"]

    if haskey(bus, "va") && !haskey(bus, "vm")
        constraint_mc_theta_ref(pm, id, nw=mw, bounded=bounded, report=report, kwargs...)
    elseif haskey(bus, "vm") && !haskey(bus, "va")
        error("Reference buses with magnitude-only (vm) setpoints are not supported. The same can be achieved by setting vmin=vmax=vm.")
    elseif haskey(bus, "vm") && haskey(bus, "va")
        constraint_mc_voltage_fixed(pm, nw, id, bus["vm"], bus["va"], terminals, grounded)
    end
end


"""
    function constraint_mc_generator_power(
        pm::ExplicitNeutralModels,
        id::Int;
        nw::Int=nw_id_default,
        report::Bool=true
    )

Constrains generator power variables for models with explicit neutrals.
"""
function constraint_mc_generator_power(pm::ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, report::Bool=true)
    generator = ref(pm, nw, :gen, id)
    bus = ref(pm, nw,:bus, generator["gen_bus"])

    configuration = generator["configuration"]

    N = length(generator["connections"])
    pmin = get(generator, "pmin", fill(-Inf, N))
    pmax = get(generator, "pmax", fill( Inf, N))
    qmin = get(generator, "qmin", fill(-Inf, N))
    qmax = get(generator, "qmax", fill( Inf, N))

    if configuration==WYE || length(pmin)==1
        constraint_mc_generator_power_wye(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax; report=report)
    else
        constraint_mc_generator_power_delta(pm, nw, id, bus["index"], generator["connections"], pmin, pmax, qmin, qmax; report=report)
    end
end


"""
    function constraint_mc_transformer_voltage(
        pm::ExplicitNeutralModels,
        i::Int;
        nw::Int=nw_id_default,
        fix_taps::Bool=true
    )

For models with explicit neutrals,
links the voltage of the from-side and to-side transformer windings
"""
function constraint_mc_transformer_voltage(pm::ExplicitNeutralModels, i::Int; nw::Int=nw_id_default, fix_taps::Bool=true)
    transformer = ref(pm, nw, :transformer, i)
    f_bus = transformer["f_bus"]
    t_bus = transformer["t_bus"]
    f_idx = (i, f_bus, t_bus)
    t_idx = (i, t_bus, f_bus)
    configuration = transformer["configuration"]
    f_connections = transformer["f_connections"]
    t_connections = transformer["t_connections"]
    tm_set = transformer["tm_set"]
    tm_fixed = fix_taps ? ones(Bool, length(tm_set)) : transformer["tm_fix"]
    tm_scale = calculate_tm_scale(transformer, ref(pm, nw, :bus, f_bus), ref(pm, nw, :bus, t_bus))

    #TODO change data model
    # there is redundancy in specifying polarity seperately on from and to side
    #TODO change this once migrated to new data model
    pol = transformer["polarity"]

    if configuration == WYE
        constraint_mc_transformer_voltage_yy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == DELTA
        constraint_mc_transformer_voltage_dy(pm, nw, i, f_bus, t_bus, f_idx, t_idx, f_connections, t_connections, pol, tm_set, tm_fixed, tm_scale)
    elseif configuration == "zig-zag"
        error("Zig-zag not yet supported.")
    end
end


"""
	function constraint_mc_transformer_thermal_limit(
		pm::ExplicitNeutralModels,
		id::Int;
		nw::Int=nw_id_default,
		bounded::Bool=true,
		report::Bool=true,
		kwargs...
	)

Imposes a bound on the total apparent at each transformer winding
"""
function constraint_mc_transformer_thermal_limit(pm::ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true, kwargs...)
    trans = ref(pm, nw, :transformer, id)
    f_bus = trans["f_bus"]
    t_bus = trans["t_bus"]
    f_idx = (id,f_bus,t_bus)
    t_idx = (id,t_bus,f_bus)
    f_conns = trans["f_connections"]
    t_conns = trans["t_connections"]
    config = trans["configuration"]
    sm_ub = trans["sm_ub"]

    constraint_mc_transformer_thermal_limit(pm, nw, id, f_idx, t_idx, f_bus, t_bus, f_conns, t_conns, config, sm_ub)
end


"""
    function constraint_mc_switch_power(
        pm::ExplicitNeutralModels,
        id::Int;
        nw::Int=nw_id_default,
        report::Bool=true
    )

For IVR models with explicit neutrals,
link the switch power or create appropiate expressions for them
"""
function constraint_mc_switch_power(pm::ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, report::Bool=true)
    switch = ref(pm, nw, :switch, id)
    f_bus = switch["f_bus"]
    t_bus = switch["t_bus"]
    f_idx = (id, f_bus, t_bus)
    t_idx = (id, t_bus, f_bus)

    constraint_mc_switch_power(pm, nw, id, f_idx, t_idx, switch["f_connections"], switch["t_connections"])
end


"""
    function constraint_mc_switch_current(
        pm::ExplicitNeutralModels,
        id::Int;
        nw::Int=nw_id_default,
        report::Bool=true
    )

For models with explicit neutrals,
link the switch currents or create appropiate expressions for them.
"""
function constraint_mc_switch_current(pm::ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, report::Bool=true)
    switch = ref(pm, nw, :switch, id)
    f_bus = switch["f_bus"]
    t_bus = switch["t_bus"]
    f_idx = (id, f_bus, t_bus)
    t_idx = (id, t_bus, f_bus)

    constraint_mc_switch_current(pm, nw, id, f_idx, t_idx, switch["f_connections"], switch["t_connections"])
end


"""
	function constraint_mc_load_power(
		pm::ExplicitNeutralModels,
		id::Int;
		nw::Int=nw_id_default,
		report::Bool=true
	)

Constrains load power variables for models with explicit neutrals.
"""
function constraint_mc_load_power(pm::ExplicitNeutralModels, id::Int; nw::Int=nw_id_default, report::Bool=true)
    load = ref(pm, nw, :load, id)
    bus = ref(pm, nw,:bus, load["load_bus"])

    configuration = load["configuration"]

    a, alpha, b, beta = _load_expmodel_params(load, bus)

    if configuration==WYE || length(a)==1
        constraint_mc_load_power_wye(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    else
        constraint_mc_load_power_delta(pm, nw, id, load["load_bus"], load["connections"], a, alpha, b, beta; report=report)
    end
end
