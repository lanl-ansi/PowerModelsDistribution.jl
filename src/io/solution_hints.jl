"""
    parse_solution(solution_file, data)

Reads a PowerModelsDistribution JSON solution and normalizes it for use as
starting values or as a JuMP feasibility point.
"""
function parse_solution(
    solution_file::String,
    data::Union{AbstractDict{String,<:Any},String},
)::Dict
    @assert endswith(lowercase(solution_file), ".json") """
    Only JSON solution files are supported
    """

    result = JSON.parsefile(
        solution_file;
        dicttype = Dict{String,Any},
    )

    return parse_solution(result, data)
end


"""
    parse_solution(result, data)

Normalizes a PowerModelsDistribution result dictionary and converts the
solution to per-unit values.
"""
function parse_solution(
    result::AbstractDict,
    data::Union{AbstractDict{String,<:Any},String},
)::Dict
    result = deepcopy(result)

    solution = normalize_solution_base_values!(result, data)

    if !get(solution, "per_unit", false)
        make_per_unit!(solution)
    end

    if haskey(result, "objective") && !haskey(solution, "objective")
        solution["objective"] = result["objective"]
    end

    return solution
end


const _PMD_SOLUTION_COMPONENTS = Set([
    "bus",
    "branch",
    "switch",
    "transformer",
    "load",
    "shunt",
    "gen",
    "storage",
])


"""
    _is_solution_root(solution)

Returns true when `solution` appears to be a PowerModelsDistribution solution
dictionary.
"""
function _is_solution_root(solution::AbstractDict)::Bool
    return haskey(solution, "nw") ||
           any(haskey(solution, component) for component in _PMD_SOLUTION_COMPONENTS)
end


"""
    normalize_solution_base_values!(result, data)

Normalizes a PowerModelsDistribution solution to SI units, replaces its base
values with those from `data`, and returns the solution root.
"""
function normalize_solution_base_values!(
    result::AbstractDict,
    data::Union{AbstractDict{String,<:Any},String},
)::Dict
    if data isa String
        data = PowerModelsDistribution.parse_file(data)
    end

    sol = _solution_root(result)

    @assert _is_solution_root(sol) "Unsupported PowerModelsDistribution solution dictionary."

    if !get(sol, "si_units", false)
        make_si_units!(sol)
    end

    @assert get(sol, "si_units", false) """
    Solution data could not be normalized to SI units.
    """

    @assert !get(sol, "per_unit", false) """
    Solution data could not be normalized to SI units.
    """

    pmd_data = get_pmd_data(data)
    _replace_base_values!(sol, pmd_data)

    return sol
end


"""
    _replace_base_values!(solution, data)

Removes base values stored in `solution` and replaces them with base values
from the supplied PowerModelsDistribution data.
"""
function _replace_base_values!(
    solution::AbstractDict,
    data::AbstractDict,
)
    for key in collect(keys(solution))
        _is_pmd_base_value(key) && delete!(solution, key)
    end

    for key in keys(data)
        _is_pmd_base_value(key) && (solution[key] = deepcopy(data[key]))
    end

    return solution
end


"""
    _is_pmd_base_value(key)

Identifies common PowerModelsDistribution base-value fields.

The `base_` check is retained for compatibility with older or custom PMD data
formats. Mathematical PMD datasets commonly use `baseMVA`.
"""
function _is_pmd_base_value(key)::Bool
    key_string = String(key)

    return startswith(key_string, "base_") ||
           key_string in ("baseMVA", "base_kv", "base_frequency")
end


"""
    add_solution_hints!(case, solution_file)

Adds values from a PowerModelsDistribution solution file as starting values in
`case`. Each solution variable `x` is added to the corresponding component as
`x_start`.
"""
function add_solution_hints!(
    case::AbstractDict,
    solution_file::String,
)::Dict
    sol_root = parse_solution(solution_file, case)

    case_root = get_pmd_data(case)

    if haskey(sol_root, "nw")
        @assert haskey(case_root, "nw") """
        The solution is multinetwork, but the destination case is not.
        """

        for (nw_id, nw_sol) in sol_root["nw"]
            case_nw_key = _lookup_template_key(case_root["nw"], nw_id)
            isnothing(case_nw_key) && continue

            _apply_hints!(nw_sol, case_root["nw"][case_nw_key])
        end
    else
        _apply_hints!(sol_root, case_root)
    end

    return case
end


"""
    build_solution_point(pm, solution)

Builds the variable-value dictionary required by
`JuMP.primal_feasibility_report`.
"""
function build_solution_point(
    pm::_IM.AbstractInfrastructureModel,
    solution::AbstractDict,
)::Dict{JuMP.VariableRef,Float64}
    point = Dict{JuMP.VariableRef,Float64}()

    objective_value = _solution_objective(solution)
    sol_root = _solution_root(solution)

    nw_sol = pm.sol[:it][pmd_it_sym][:nw]

    if haskey(sol_root, "nw")
        for (nw_id, nw_src) in sol_root["nw"]
            nw_key = _lookup_template_key(nw_sol, nw_id)
            isnothing(nw_key) && continue

            _build_solution_point!(
                point,
                nw_sol[nw_key],
                nw_src,
            )
        end
    else
        nw_key = _lookup_template_key(nw_sol, nw_id_default)

        if isnothing(nw_key)
            error(
                "Could not find default network $(nw_id_default) " *
                "in the PowerModelsDistribution solution template.",
            )
        end

        _build_solution_point!(
            point,
            nw_sol[nw_key],
            sol_root,
        )
    end

    _add_objective_aux_value!(
        point,
        pm.model,
        objective_value,
    )

    return point
end


"""
    primal_feasibility_report(pm, solution; atol=0.0, skip_missing=false)

Creates a JuMP primal-feasibility report from a PMD result dictionary.
"""
function primal_feasibility_report(
    pm::_IM.AbstractInfrastructureModel,
    solution::AbstractDict;
    atol::Float64 = 0.0,
    skip_missing::Bool = false,
)
    normalized_solution = parse_solution(solution, pm.data)
    point = build_solution_point(pm, normalized_solution)

    return JuMP.primal_feasibility_report(
        pm.model,
        point;
        atol = atol,
        skip_missing = skip_missing,
    )
end


"""
    primal_feasibility_report(pm, solution_file; atol=0.0, skip_missing=false)

Creates a JuMP primal-feasibility report from a PMD JSON solution file.
"""
function primal_feasibility_report(
    pm::_IM.AbstractInfrastructureModel,
    solution_file::String;
    atol::Float64 = 0.0,
    skip_missing::Bool = false,
)
    normalized_solution = parse_solution(solution_file, pm.data)
    point = build_solution_point(pm, normalized_solution)

    return JuMP.primal_feasibility_report(
        pm.model,
        point;
        atol = atol,
        skip_missing = skip_missing,
    )
end


"""
    _apply_hints!(src, dst)

Copies scalar solution values into corresponding `_start` fields in a
PowerModelsDistribution data dictionary.
"""
function _apply_hints!(
    src::AbstractDict,
    dst::AbstractDict,
)
    for (component_name, component_solution) in src
        component_solution isa AbstractDict || continue

        dst_component_key = _lookup_template_key(dst, component_name)
        isnothing(dst_component_key) && continue

        destination_components = dst[dst_component_key]
        destination_components isa AbstractDict || continue

        for (component_id, component_values) in component_solution
            component_values isa AbstractDict || continue

            dst_id_key = _lookup_template_key(
                destination_components,
                component_id,
            )
            isnothing(dst_id_key) && continue

            destination_component = destination_components[dst_id_key]
            destination_component isa AbstractDict || continue

            for (variable_name, value) in component_values
                value isa Number || continue

                start_key = string(variable_name, "_start")
                destination_component[start_key] = Float64(value)
            end
        end
    end

    return dst
end


"""
    _solution_root(solution)

Returns the inner `"solution"` dictionary when present.
"""
function _solution_root(solution::AbstractDict)::AbstractDict
    if haskey(solution, "solution")
        solution["solution"] isa AbstractDict || error(
            "Unsupported solution format: `solution` must map to a dictionary",
        )

        return solution["solution"]
    end

    return solution
end


"""
    _solution_objective(solution)

Returns an objective value stored either on the outer result dictionary or
inside its solution root.
"""
function _solution_objective(solution::AbstractDict)
    if haskey(solution, "objective") && solution["objective"] isa Number
        return solution["objective"]
    end

    sol_root = _solution_root(solution)

    if haskey(sol_root, "objective") && sol_root["objective"] isa Number
        return sol_root["objective"]
    end

    return nothing
end


"""
    _build_solution_point!(point, template, value)

Recursively matches values in a serialized PMD solution to the corresponding
JuMP variables in the PMD solution template.
"""
function _build_solution_point!(
    point::AbstractDict{JuMP.VariableRef,Float64},
    template,
    value,
)
    if template isa JuMP.VariableRef && value isa Number
        point[template] = Float64(value)

    elseif template isa AbstractArray && value isa AbstractArray
        for index in eachindex(template, value)
            _build_solution_point!(
                point,
                template[index],
                value[index],
            )
        end

    elseif template isa AbstractDict && value isa AbstractDict
        for (key, subvalue) in value
            template_key = _lookup_template_key(template, key)
            isnothing(template_key) && continue

            _build_solution_point!(
                point,
                template[template_key],
                subvalue,
            )
        end
    end

    return point
end


"""
    _lookup_template_key(template, key)

Matches string, integer, and symbol representations of dictionary keys.
"""
function _lookup_template_key(
    template::AbstractDict,
    key,
)
    if haskey(template, key)
        return key
    end

    if key isa AbstractString
        symbol_key = Symbol(key)

        if haskey(template, symbol_key)
            return symbol_key
        end

        integer_key = tryparse(Int, key)

        if !isnothing(integer_key) && haskey(template, integer_key)
            return integer_key
        end
    elseif key isa Symbol
        string_key = String(key)

        if haskey(template, string_key)
            return string_key
        end
    elseif key isa Integer
        string_key = string(key)

        if haskey(template, string_key)
            return string_key
        end
    end

    return nothing
end


"""
    _add_objective_aux_value!(point, model, objective_value)

Adds the supplied objective value when the JuMP objective is represented by a
single auxiliary variable.
"""
function _add_objective_aux_value!(
    point::AbstractDict{JuMP.VariableRef,Float64},
    model::JuMP.Model,
    objective_value,
)
    objective_value isa Number || return point

    objective_function = JuMP.objective_function(model)

    if objective_function isa JuMP.VariableRef &&
       !haskey(point, objective_function)
        point[objective_function] = Float64(objective_value)
    end

    return point
end