"on-load tap-changer OPF with ACPPowerModel"
function run_ac_mc_opf_oltc(data::Union{Dict{String,<:Any},String}, solver; kwargs...)
    return run_mc_opf_oltc(data, ACPPowerModel, solver; kwargs...)
end


"on-load tap-changer OPF"
function run_mc_opf_oltc(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_opf_oltc; kwargs...)
end


"constructor for on-load tap-changer OPF"
function build_mc_opf_oltc(pm::_PM.AbstractPowerModel)
    variable_mc_bus_voltage(pm)
    variable_mc_branch_power(pm)
    variable_mc_gen_power_setpoint(pm)
    variable_mc_load_setpoint(pm)
    variable_mc_transformer_power(pm)
    variable_mc_oltc_transformer_tap(pm)

    constraint_mc_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    # generators should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_gen_setpoint(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_setpoint(pm, id)
    end

    for i in ids(pm, :bus)
        constraint_mc_load_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i, fix_taps=false)
    end

    _PM.objective_min_fuel_cost(pm)
end
