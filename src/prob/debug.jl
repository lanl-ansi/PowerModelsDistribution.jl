# These problem formulations are used to debug Distribution datasets
# that do not converge using the standard formulations
"Solve OPF problem with slack power at every bus"
function solve_mc_opf_pbs(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_opf_pbs; kwargs...)
end


"Solve PF problem with slack power at every bus"
function solve_mc_pf_pbs(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_pf_pbs; kwargs...)
end


"OPF problem with slack power at every bus"
function build_mc_opf_pbs(pm::_PM.AbstractPowerModel)
    variable_mc_bus_voltage(pm)

    variable_mc_branch_power(pm)
    variable_mc_switch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_generator_power(pm)

    variable_mc_slack_bus_power(pm)

    constraint_mc_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        constraint_mc_power_balance_slack(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    objective_mc_min_slack_bus_power(pm)
end


"PF problem with slack power at every bus"
function build_mc_pf_pbs(pm::_PM.AbstractPowerModel)
    variable_mc_bus_voltage(pm; bounded=false)

    variable_mc_branch_power(pm; bounded=false)
    variable_mc_switch_power(pm; bounded=false)
    variable_mc_transformer_power(pm; bounded=false)

    variable_mc_generator_power(pm; bounded=false)

    variable_mc_slack_bus_power(pm)

    constraint_mc_model_voltage(pm)

    for (i,bus) in ref(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)

        @assert bus["bus_type"] == 3
        constraint_mc_voltage_magnitude_only(pm, i)
    end

    for (i,bus) in ref(pm, :bus)
        constraint_mc_power_balance_slack(pm, i)

        # PV Bus Constraints
        if length(ref(pm, :bus_gens, i)) > 0 && !(i in ids(pm,:ref_buses))
            # this assumes inactive generators are filtered out of bus_gens
            @assert bus["bus_type"] == 2

            constraint_mc_voltage_magnitude_only(pm, i)
            for j in ref(pm, :bus_gens, i)
                constraint_mc_gen_power_setpoint_real(pm, j)
            end
        end
    end

    for i in ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)
    end

    objective_mc_min_slack_bus_power(pm)
end

# Depreciated run_ functions (remove after ~4-6 months)

"depreciation message for run_mc_opf_pbs"
function run_mc_opf_pbs(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    @warn "run_mc_opf_pbs is being depreciated in favor of solve_mc_opf_pbs, please update your code accordingly"
    return solve_mc_opf_pbs(data, model_type, solver; kwargs...)
end


"depreciation message for run_mc_opf_pbs"
function run_mc_pf_pbs(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    @warn "run_mc_pf_pbs is being depreciated in favor of solve_mc_pf_pbs, please update your code accordingly"
    return solve_mc_pf_pbs(data, model_type, solver; kwargs...)
end
