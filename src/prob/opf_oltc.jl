"Solve on-load tap-changer OPF"
function solve_mc_opf_oltc(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mc_opf_oltc; kwargs...)
end


"constructor for on-load tap-changer OPF"
function build_mc_opf_oltc(pm::AbstractUnbalancedPowerModel)
    variable_mc_bus_voltage(pm)

    variable_mc_branch_power(pm)
    variable_mc_switch_power(pm)
    variable_mc_transformer_power(pm)

    variable_mc_oltc_transformer_tap(pm)

    variable_mc_generator_power(pm)
    variable_mc_load_power(pm)

    constraint_mc_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        constraint_mc_theta_ref(pm, i)
    end

    # generators should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :gen)
        constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in ids(pm, :load)
        constraint_mc_load_power(pm, id)
    end

    for i in ids(pm, :bus)
        constraint_mc_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        constraint_mc_ohms_yt_from(pm, i)
        constraint_mc_ohms_yt_to(pm, i)

        constraint_mc_voltage_angle_difference(pm, i)

        constraint_mc_thermal_limit_from(pm, i)
        constraint_mc_thermal_limit_to(pm, i)
    end

    for i in ids(pm, :switch)
        constraint_mc_switch_state(pm, i)
        constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in ids(pm, :transformer)
        constraint_mc_transformer_power(pm, i, fix_taps=false)
    end

    objective_mc_min_fuel_cost(pm)
end

# Depreciated run_ functions (remove after ~4-6 months)

"depreciation warning for run_ac_mc_opf_oltc"
function run_ac_mc_opf_oltc(data::Union{Dict{String,<:Any},String}, solver; kwargs...)
    @warn "run_ac_mc_opf_oltc is being depreciated in favor of solve_mc_opf_oltc(data, ACPUPowerModel, solver; kwargs...), please update your code accordingly"
    return solve_mc_opf_oltc(data, ACPUPowerModel, solver; kwargs...)
end


"depreciation warning for run_mc_opf_oltc"
function run_mc_opf_oltc(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    @warn "run_mc_opf_oltc is being depreciated in favor of solve_mc_opf_oltc, please update your code accordingly"
    return solve_mc_opf_oltc(data, model_type, solver; kwargs...)
end
