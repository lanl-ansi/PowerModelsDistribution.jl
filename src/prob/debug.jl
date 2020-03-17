# These problem formulations are used to debug Distribution datasets
# that do not converge using the standard formulations
"OPF problem with slack power at every bus"
function run_mc_opf_pbs(data::Dict{String,Any}, model_type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_opf_pbs; kwargs...)
end


"OPF problem with slack power at every bus"
function run_mc_opf_pbs(file::String, model_type, solver; kwargs...)
    return run_mc_opf_pbs(PowerModelsDistribution.parse_file(file), model_type, solver; kwargs...)
end


"PF problem with slack power at every bus"
function run_mc_pf_pbs(data::Dict{String,Any}, model_type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, build_mc_pf_pbs; kwargs...)
end


"PF problem with slack power at every bus"
function run_mc_pf_pbs(file::String, model_type, solver; kwargs...)
    return run_mc_pf_pbs(PowerModelsDistribution.parse_file(file), model_type, solver; kwargs...)
end


"OPF problem with slack power at every bus"
function build_mc_opf_pbs(pm::_PMs.AbstractPowerModel)
    variable_mc_voltage(pm)

    variable_mc_branch_flow(pm)
    variable_mc_transformer_flow(pm)

    variable_mc_generation(pm)

    variable_mc_bus_power_slack(pm)

    constraint_mc_model_voltage(pm)

    for i in _PMs.ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    for i in _PMs.ids(pm, :bus)
        constraint_mc_power_balance_slack(pm, i)
    end

    for i in _PMs.ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    objective_min_bus_power_slack(pm)
end


"PF problem with slack power at every bus"
function build_mc_pf_pbs(pm::_PMs.AbstractPowerModel)
    variable_mc_voltage(pm; bounded=false)

    variable_mc_branch_flow(pm; bounded=false)
    variable_mc_transformer_flow(pm; bounded=false)

    variable_mc_generation(pm; bounded=false)

    variable_mc_bus_power_slack(pm)

    constraint_mc_model_voltage(pm)

    for (i,bus) in _PMs.ref(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)

        @assert bus["bus_type"] == 3
        constraint_mc_voltage_magnitude_setpoint(pm, i)
    end

    for (i,bus) in _PMs.ref(pm, :bus)
        constraint_mc_power_balance_slack(pm, i)

        # PV Bus Constraints
        if length(_PMs.ref(pm, :bus_gens, i)) > 0 && !(i in _PMs.ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_mc_voltage_magnitude_setpoint(pm, i)
            for j in _PMs.ref(pm, :bus_gens, i)
                constraint_mc_active_gen_setpoint(pm, j)
            end
        end
    end

    for i in _PMs.ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)
    end

    objective_min_bus_power_slack(pm)
end
