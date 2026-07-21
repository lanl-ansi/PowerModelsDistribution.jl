"""
    comprehensive_feasibility_report(
        data::MathematicalModel,
        model_type::Type;
        build_method = build_mc_opf,
        optimizer = nothing,
        point_source::Symbol = :start,
        solve::Bool = false,
        atol::Float64 = 1e-6,
        top::Int = 25,
        fill_missing::Symbol = :error,
        bound_atol::Float64 = atol,
        large_value_threshold::Float64 = 1e4,
        io::IO = stdout,
        instantiate_kwargs...,
    )

Instantiate a PowerModelsDistribution optimization model and generate a
detailed primal-feasibility report.

# Point sources

- `point_source = :start` evaluates the JuMP variable start values.
- `point_source = :solution` evaluates the most recent solver result.

Set `solve = true` to optimize the model before generating the report.

# Missing start handling

When `point_source = :start`, `fill_missing` may be:

- `:error`: throw an error if any variables lack starts.
- `:zero`: use zero for missing starts, projected onto variable bounds.
- `:bounds`: use a fixed value, bound midpoint, or nearest finite bound.
- `:skip`: omit missing variables and ask JuMP to skip constraints that use
  them. This is incompatible with some legacy nonlinear constraints.

# Returns

A named tuple containing:

- `pm`: instantiated PowerModelsDistribution model.
- `solve_result`: result from `optimize_model!`, or `nothing`.
- `point`: variable-value dictionary evaluated by JuMP.
- `raw_report`: raw output from `JuMP.primal_feasibility_report`.
- `summary`: aggregate diagnostic information.
- `violations`: sorted constraint-violation records.
"""
function comprehensive_feasibility_report(
    data::MathematicalModel,
    model_type::Type;
    build_method = build_mc_opf,
    optimizer = nothing,
    point_source::Symbol = :start,
    solve::Bool = false,
    atol::Float64 = 1e-6,
    top::Int = 10,
    fill_missing::Symbol = :error,
    bound_atol::Float64 = atol,
    large_value_threshold::Float64 = 1e4,
    io::IO = stdout,
    instantiate_kwargs...,
)
    point_source in (:start, :solution) || throw(
        ArgumentError(
            "`point_source` must be either `:start` or `:solution`, " *
            "not `$(point_source)`.",
        ),
    )

    fill_missing in (:error, :zero, :bounds, :skip) || throw(
        ArgumentError(
            "`fill_missing` must be one of `:error`, `:zero`, " *
            "`:bounds`, or `:skip`.",
        ),
    )

    top >= 0 || throw(ArgumentError("`top` must be nonnegative."))
    atol >= 0.0 || throw(ArgumentError("`atol` must be nonnegative."))
    bound_atol >= 0.0 || throw(
        ArgumentError("`bound_atol` must be nonnegative."),
    )

    pm = instantiate_mc_model(
        data,
        model_type,
        build_method;
        instantiate_kwargs...,
    )

    solve_result = nothing

    if solve
        isnothing(optimizer) && throw(
            ArgumentError(
                "An optimizer must be supplied when `solve = true`.",
            ),
        )

        solve_result = optimize_model!(
            pm;
            optimizer = optimizer,
        )
    elseif !isnothing(optimizer)
        JuMP.set_optimizer(pm.model, optimizer)
    end

    model = pm.model
    variables = JuMP.all_variables(model)

    point, point_diagnostics = _build_diagnostic_point(
        variables,
        point_source;
        fill_missing = fill_missing,
    )

    skip_missing = fill_missing == :skip

    raw_report = JuMP.primal_feasibility_report(
        model,
        point;
        atol = atol,
        skip_missing = skip_missing,
    )

    violations = _build_violation_records(
        raw_report,
        point,
    )

    variable_summary = _summarize_variables(
        variables,
        point;
        bound_atol = bound_atol,
        large_value_threshold = large_value_threshold,
    )

    constraint_summary = _summarize_constraints(
        model,
        violations,
    )

    solver_summary = _solver_summary(
        model,
        point_source,
    )

    summary = (
        point_source = point_source,
        tolerance = atol,
        variable = variable_summary,
        constraint = constraint_summary,
        solver = solver_summary,
        point = point_diagnostics,
    )

    _print_feasibility_report(
        io,
        model_type,
        summary,
        violations;
        top = top,
    )

    return (
        pm = pm,
        solve_result = solve_result,
        point = point,
        raw_report = raw_report,
        summary = summary,
        violations = violations,
    )
end


"""
Build the point supplied to `JuMP.primal_feasibility_report`.
"""
function _build_diagnostic_point(
    variables::Vector{<:JuMP.VariableRef},
    point_source::Symbol;
    fill_missing::Symbol,
)
    point = Dict{JuMP.VariableRef,Float64}()

    missing_variables = JuMP.VariableRef[]
    nonfinite_variables = Pair{JuMP.VariableRef,Float64}[]
    filled_variables = Pair{JuMP.VariableRef,Float64}[]

    for variable in variables
        raw_value = if point_source == :start
            JuMP.start_value(variable)
        else
            JuMP.has_values(JuMP.owner_model(variable)) || error(
                "The model has no primal solution. Set `solve = true` or " *
                "use `point_source = :start`.",
            )

            JuMP.value(variable)
        end

        if isnothing(raw_value)
            push!(missing_variables, variable)

            if fill_missing == :error || fill_missing == :skip
                continue
            end

            filled_value = _default_variable_value(
                variable,
                fill_missing,
            )

            point[variable] = filled_value
            push!(filled_variables, variable => filled_value)
            continue
        end

        value = Float64(raw_value)

        if !isfinite(value)
            push!(nonfinite_variables, variable => value)
            continue
        end

        point[variable] = value
    end

    if !isempty(missing_variables) && fill_missing == :error
        variable_names = join(
            ["  - " * _variable_label(v) for v in missing_variables],
            "\n",
        )

        error(
            "$(length(missing_variables)) variables have no start value:\n" *
            variable_names *
            "\nUse `fill_missing = :zero`, `:bounds`, or `:skip` " *
            "to continue.",
        )
    end

    if !isempty(nonfinite_variables)
        variable_names = join(
            [
                "  - $(_variable_label(v)) = $(value)"
                for (v, value) in nonfinite_variables
            ],
            "\n",
        )

        error(
            "$(length(nonfinite_variables)) variables have nonfinite " *
            "values:\n$(variable_names)",
        )
    end

    diagnostics = (
        missing = missing_variables,
        filled = filled_variables,
        nonfinite = nonfinite_variables,
    )

    return point, diagnostics
end


"""
Choose a finite replacement for a missing variable start.
"""
function _default_variable_value(
    variable::JuMP.VariableRef,
    strategy::Symbol,
)::Float64
    if JuMP.is_fixed(variable)
        return Float64(JuMP.fix_value(variable))
    end

    has_lower = JuMP.has_lower_bound(variable)
    has_upper = JuMP.has_upper_bound(variable)

    lower = has_lower ? Float64(JuMP.lower_bound(variable)) : -Inf
    upper = has_upper ? Float64(JuMP.upper_bound(variable)) : Inf

    value = if strategy == :zero
        0.0
    elseif strategy == :bounds
        if isfinite(lower) && isfinite(upper)
            (lower + upper) / 2.0
        elseif isfinite(lower)
            lower
        elseif isfinite(upper)
            upper
        else
            0.0
        end
    else
        error("Unsupported missing-value strategy: $(strategy)")
    end

    return clamp(value, lower, upper)
end


"""
Convert JuMP's raw violation dictionary into sorted diagnostic records.
"""
function _build_violation_records(
    raw_report::AbstractDict,
    point::AbstractDict{JuMP.VariableRef,Float64},
)
    records = NamedTuple[]

    for (constraint, violation) in raw_report
        object = JuMP.constraint_object(constraint)

        function_value = try
            JuMP.value(
                variable -> get(point, variable, missing),
                object.func,
            )
        catch
            missing
        end

        constraint_name = try
            JuMP.name(constraint)
        catch
            ""
        end

        push!(
            records,
            (
                constraint = constraint,
                violation = Float64(violation),
                name = constraint_name,
                label = _constraint_label(constraint),
                function_type = typeof(object.func),
                set_type = typeof(object.set),
                function_value = function_value,
                set = object.set,
                is_variable_bound = object.func isa JuMP.VariableRef,
            ),
        )
    end

    sort!(
        records;
        by = record -> record.violation,
        rev = true,
    )

    return records
end


"""
Create aggregate variable diagnostics.
"""
function _summarize_variables(
    variables::Vector{<:JuMP.VariableRef},
    point::AbstractDict{JuMP.VariableRef,Float64};
    bound_atol::Float64,
    large_value_threshold::Float64,
)
    point_values = collect(Base.values(point))

    fixed = JuMP.VariableRef[]
    at_lower_bound = Pair{JuMP.VariableRef,Float64}[]
    at_upper_bound = Pair{JuMP.VariableRef,Float64}[]
    lower_bound_violations = NamedTuple[]
    upper_bound_violations = NamedTuple[]
    fixed_value_violations = NamedTuple[]
    zero_values = JuMP.VariableRef[]
    large_values = Pair{JuMP.VariableRef,Float64}[]

    for variable in variables
        haskey(point, variable) || continue

        value = point[variable]

        if JuMP.is_fixed(variable)
            push!(fixed, variable)

            fixed_value = Float64(JuMP.fix_value(variable))
            violation = abs(value - fixed_value)

            if violation > bound_atol
                push!(
                    fixed_value_violations,
                    (
                        variable = variable,
                        value = value,
                        fixed_value = fixed_value,
                        violation = violation,
                    ),
                )
            end
        end

        if JuMP.has_lower_bound(variable)
            lower = Float64(JuMP.lower_bound(variable))
            distance = value - lower

            if value < lower - bound_atol
                push!(
                    lower_bound_violations,
                    (
                        variable = variable,
                        value = value,
                        bound = lower,
                        violation = lower - value,
                    ),
                )
            elseif abs(distance) <= bound_atol
                push!(at_lower_bound, variable => value)
            end
        end

        if JuMP.has_upper_bound(variable)
            upper = Float64(JuMP.upper_bound(variable))
            distance = upper - value

            if value > upper + bound_atol
                push!(
                    upper_bound_violations,
                    (
                        variable = variable,
                        value = value,
                        bound = upper,
                        violation = value - upper,
                    ),
                )
            elseif abs(distance) <= bound_atol
                push!(at_upper_bound, variable => value)
            end
        end

        iszero(value) && push!(zero_values, variable)

        if abs(value) >= large_value_threshold
            push!(large_values, variable => value)
        end
    end

    sort!(
        lower_bound_violations;
        by = item -> item.violation,
        rev = true,
    )

    sort!(
        upper_bound_violations;
        by = item -> item.violation,
        rev = true,
    )

    sort!(
        fixed_value_violations;
        by = item -> item.violation,
        rev = true,
    )

    sort!(
        large_values;
        by = item -> abs(last(item)),
        rev = true,
    )

    return (
        total = length(variables),
        represented_in_point = length(point),
        fixed = fixed,
        zero_values = zero_values,
        at_lower_bound = at_lower_bound,
        at_upper_bound = at_upper_bound,
        lower_bound_violations = lower_bound_violations,
        upper_bound_violations = upper_bound_violations,
        fixed_value_violations = fixed_value_violations,
        large_values = large_values,
        minimum = isempty(point_values) ? missing : minimum(point_values),
        maximum = isempty(point_values) ? missing : maximum(point_values),
        maximum_absolute = isempty(point_values) ?
            missing :
            maximum(abs, point_values),
        two_norm = isempty(point_values) ?
            missing :
            sqrt(sum(abs2, point_values)),
    )
end


"""
Summarize all model constraints and the violated subset by type.
"""
function _summarize_constraints(
    model::JuMP.Model,
    violations::Vector{<:NamedTuple},
)
    total_by_type = Dict{Tuple{DataType,DataType},Int}()

    total_constraints = 0

    for (function_type, set_type) in JuMP.list_of_constraint_types(model)
        count = JuMP.num_constraints(
            model,
            function_type,
            set_type,
        )

        total_by_type[(function_type, set_type)] = count
        total_constraints += count
    end

    violated_by_type = Dict{Tuple{DataType,DataType},Int}()

    for record in violations
        key = (record.function_type, record.set_type)
        violated_by_type[key] = get(violated_by_type, key, 0) + 1
    end

    violation_values = [record.violation for record in violations]

    return (
        total = total_constraints,
        violated = length(violations),
        satisfied = total_constraints - length(violations),
        total_by_type = total_by_type,
        violated_by_type = violated_by_type,
        maximum_violation = isempty(violation_values) ?
            0.0 :
            maximum(violation_values),
        mean_violation = isempty(violation_values) ?
            0.0 :
            sum(violation_values) / length(violation_values),
        violation_two_norm = sqrt(sum(abs2, violation_values)),
        bound_violations = count(
            record -> record.is_variable_bound,
            violations,
        ),
        nonbound_violations = count(
            record -> !record.is_variable_bound,
            violations,
        ),
    )
end


"""
Collect solver and objective information without assuming the model solved.
"""
function _solver_summary(
    model::JuMP.Model,
    point_source::Symbol,
)
    has_solution = JuMP.has_values(model)

    termination_status = has_solution ?
        JuMP.termination_status(model) :
        missing

    primal_status = has_solution ?
        JuMP.primal_status(model) :
        missing

    dual_status = has_solution ?
        JuMP.dual_status(model) :
        missing

    objective_value = if has_solution
        try
            JuMP.objective_value(model)
        catch
            missing
        end
    else
        missing
    end

    objective_at_point = if point_source == :start
        try
            JuMP.value(
                variable -> JuMP.start_value(variable),
                JuMP.objective_function(model),
            )
        catch
            missing
        end
    else
        objective_value
    end

    return (
        has_solution = has_solution,
        termination_status = termination_status,
        primal_status = primal_status,
        dual_status = dual_status,
        objective_value = objective_value,
        objective_at_point = objective_at_point,
    )
end


"""
Print the human-readable feasibility report.
"""
function _print_feasibility_report(
    io::IO,
    model_type::Type,
    summary::NamedTuple,
    violations::Vector{<:NamedTuple};
    top::Int,
)
    variable = summary.variable
    constraint = summary.constraint
    solver = summary.solver
    point = summary.point

    println(io)
    println(io, repeat("=", 80))
    println(io, "COMPREHENSIVE PRIMAL FEASIBILITY REPORT")
    println(io, repeat("=", 80))
    println(io, "Formulation:       ", model_type)
    println(io, "Point source:      ", summary.point_source)
    println(io, "Tolerance:         ", summary.tolerance)

    if !ismissing(solver.termination_status)
        println(io, "Termination:       ", solver.termination_status)
        println(io, "Primal status:     ", solver.primal_status)
        println(io, "Dual status:       ", solver.dual_status)
    end

    if !ismissing(solver.objective_at_point)
        println(io, "Objective:         ", solver.objective_at_point)
    end

    println(io)
    println(io, "VARIABLE SUMMARY")
    println(io, repeat("-", 80))
    println(io, "Total variables:           ", variable.total)
    println(io, "Variables in point:        ", variable.represented_in_point)
    println(io, "Missing values:            ", length(point.missing))
    println(io, "Filled missing values:     ", length(point.filled))
    println(io, "Fixed variables:           ", length(variable.fixed))
    println(io, "Zero-valued variables:     ", length(variable.zero_values))
    println(io, "At lower bound:            ", length(variable.at_lower_bound))
    println(io, "At upper bound:            ", length(variable.at_upper_bound))
    println(
        io,
        "Lower-bound violations:   ",
        length(variable.lower_bound_violations),
    )
    println(
        io,
        "Upper-bound violations:   ",
        length(variable.upper_bound_violations),
    )
    println(
        io,
        "Fixed-value violations:   ",
        length(variable.fixed_value_violations),
    )
    println(io, "Minimum value:             ", variable.minimum)
    println(io, "Maximum value:             ", variable.maximum)
    println(io, "Maximum absolute value:    ", variable.maximum_absolute)
    println(io, "Point two-norm:            ", variable.two_norm)

    println(io)
    println(io, "CONSTRAINT SUMMARY")
    println(io, repeat("-", 80))
    println(io, "Total constraints:         ", constraint.total)
    println(io, "Satisfied within atol:     ", constraint.satisfied)
    println(io, "Violated constraints:      ", constraint.violated)
    println(io, "Bound violations:          ", constraint.bound_violations)
    println(io, "Non-bound violations:      ", constraint.nonbound_violations)
    println(io, "Maximum violation:         ", constraint.maximum_violation)
    println(io, "Mean violated residual:    ", constraint.mean_violation)
    println(io, "Violation two-norm:        ", constraint.violation_two_norm)

    if !isempty(constraint.violated_by_type)
        println(io)
        println(io, "VIOLATIONS BY CONSTRAINT TYPE")
        println(io, repeat("-", 80))

        type_counts = collect(constraint.violated_by_type)

        sort!(
            type_counts;
            by = pair -> last(pair),
            rev = true,
        )

        for ((function_type, set_type), violated_count) in type_counts
            total_count = get(
                constraint.total_by_type,
                (function_type, set_type),
                0,
            )

            println(
                io,
                lpad(violated_count, 6),
                " / ",
                lpad(total_count, 6),
                "  ",
                function_type,
                " in ",
                set_type,
            )
        end
    end

    display_count = min(top, length(violations))

    println(io)
    println(io, "TOP $(display_count) CONSTRAINT VIOLATIONS")
    println(io, repeat("-", 80))

    if display_count == 0
        println(io, "No violations exceed the requested tolerance.")
    else
        for index in 1:display_count
            record = violations[index]

            println(
                io,
                "[",
                lpad(index, ndigits(display_count)),
                "] violation = ",
                record.violation,
            )

            println(io, "    constraint: ", record.label)

            if !isempty(record.name)
                println(io, "    name:       ", record.name)
            end

            println(io, "    function:   ", record.function_type)
            println(io, "    set:        ", record.set)

            if !ismissing(record.function_value)
                println(io, "    value:      ", record.function_value)
            end
        end
    end

    if !isempty(point.filled)
        println(io)
        println(io, "MISSING STARTS FILLED AUTOMATICALLY")
        println(io, repeat("-", 80))

        for (variable, value) in point.filled
            println(
                io,
                "  ",
                _variable_label(variable),
                " = ",
                value,
            )
        end
    elseif !isempty(point.missing)
        println(io)
        println(io, "VARIABLES MISSING VALUES")
        println(io, repeat("-", 80))

        for variable in point.missing
            println(io, "  ", _variable_label(variable))
        end
    end

    if !isempty(variable.fixed_value_violations)
        println(io)
        println(io, "FIXED-VARIABLE INCONSISTENCIES")
        println(io, repeat("-", 80))

        for item in variable.fixed_value_violations
            println(
                io,
                "  ",
                _variable_label(item.variable),
                ": value = ",
                item.value,
                ", fixed = ",
                item.fixed_value,
                ", violation = ",
                item.violation,
            )
        end
    end

    println(io, repeat("=", 80))
    println(io)

    return nothing
end


function _variable_label(variable::JuMP.VariableRef)::String
    variable_name = JuMP.name(variable)

    return isempty(variable_name) ?
        string(variable) :
        variable_name
end


function _constraint_label(constraint)::String
    constraint_name = try
        JuMP.name(constraint)
    catch
        ""
    end

    return isempty(constraint_name) ?
        string(constraint) :
        constraint_name
end