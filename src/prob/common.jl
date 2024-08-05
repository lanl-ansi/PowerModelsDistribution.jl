"""
    _solve_mc_model(
        data::Dict{String,<:Any},
        model_type::Type,
        optimizer,
        build_method::Function;
        multinetwork::Bool=false,
        ref_extensions::Vector{<:Function}=Function[],
        solution_processors::Vector{<:Function}=Function[],
        relax_integrality::Bool=false,
        kwargs...
    )::Dict{String,Any}

Internal solver interface that uses [`instantiate_mc_model`](@ref instantiate_mc_model) directly and runs `optimize_model!`, returning a result

See [`solve_mc_model`](@ref solve_mc_model)
"""
function _solve_mc_model(
    data::Dict{String,<:Any},
    model_type::Type,
    optimizer,
    build_method::Function;
    multinetwork::Bool=false,
    ref_extensions::Vector{<:Function}=Function[],
    solution_processors::Vector{<:Function}=Function[],
    relax_integrality::Bool=false,
    kwargs...
    )::Dict{String,Any}

    if multinetwork != ismultinetwork(data)
        model_requirement = multinetwork ? "multi-network" : "single-network"
        data_type = ismultinetwork(data) ? "multi-network" : "single-network"
        error("attempted to build a $(model_requirement) model with $(data_type) data")
    end

    start_time = time()
    pm = instantiate_mc_model(
        data,
        model_type,
        build_method;
        ref_extensions=ref_extensions,
        multinetwork=multinetwork,
        kwargs...
    )
    @debug "pm model build time: $(time() - start_time)"

    start_time = time()
    result = optimize_model!(
        pm,
        relax_integrality=relax_integrality,
        optimizer=optimizer,
        solution_processors=solution_processors
    )
    @debug "pm model solve and solution time: $(time() - start_time)"

    return result
end


"""
    instantiate_mc_model(
        data::Dict{String,<:Any},
        model_type::Type,
        build_method::Function;
        ref_extensions::Vector{<:Function}=Function[],
        multinetwork::Bool=false,
        global_keys::Set{String}=Set{String}(),
        eng2math_extensions::Vector{<:Function}=Function[],
        eng2math_passthrough::Dict{String,<:Vector{<:String}}=Dict{String,Vector{String}}(),
        make_pu_extensions::Vector{<:Function}=Function[],
        kwargs...
    )

Takes data in either the ENGINEERING or MATHEMATICAL model, a model type (_e.g._, [`ACRUPowerModel`](@ref ACRUPowerModel)),
and model builder function (_e.g._, [`build_mc_opf`](@ref build_mc_opf)), and returns an
[`AbstractUnbalancedPowerModel`](@ref AbstractUnbalancedPowerModel) structure.

For an explanation of `multinetwork` and `global_keys`, see [`make_multinetwork`](@ref make_multinetwork)

For an explanation of `eng2math_extensions` and `eng2math_passthrough`, see [`transform_data_model`](@ref transform_data_model)

For an explanation of `make_pu_extensions`, see [`make_per_unit!`](@ref make_per_unit!)

# `ref_extensions`

Inside of the model structures, data can be quickly accessed via the [`ref`](@ref ref) function. By default, the only ref
objects available are created by [`ref_add_core!`](@ref ref_add_core!), but users can add their own custom ref objects by passing
functions via `ref_extensions` that have the signature:

    ref_add!(ref::Dict{Symbol,Any}, data::Dict{String,Any})

See the [Beginners Guide](@ref Introduction-to-PowerModelsDistribution) for an example.
"""
function instantiate_mc_model(
    data::Dict{String,<:Any},
    model_type::Type,
    build_method::Function;
    ref_extensions::Vector{<:Function}=Function[],
    multinetwork::Bool=ismultinetwork(data),
    global_keys::Set{String}=Set{String}(),
    eng2math_extensions::Vector{<:Function}=Function[],
    eng2math_passthrough::Dict{String,<:Vector{<:String}}=Dict{String,Vector{String}}(),
    make_pu_extensions::Vector{<:Function}=Function[],
    kwargs...
    )

    if iseng(data)
        @info "Converting ENGINEERING data model to MATHEMATICAL first to build JuMP model"
        data = transform_data_model(
            data;
            multinetwork=multinetwork,
            global_keys=global_keys,
            eng2math_extensions=eng2math_extensions,
            eng2math_passthrough=eng2math_passthrough,
            make_pu_extensions=make_pu_extensions,
        )
    end

    return _IM.instantiate_model(
        data,
        model_type,
        build_method,
        ref_add_core!,
        union(_pmd_math_global_keys, global_keys),
        pmd_it_sym;
        ref_extensions=ref_extensions,
        kwargs...
    )
end



function instantiate_mc_model_ravens(
    data::Dict{String,<:Any},
    model_type::Type,
    build_method::Function;
    ref_extensions::Vector{<:Function}=Function[],
    multinetwork::Bool=ismultinetwork(data),
    global_keys::Set{String}=Set{String}(),
    ravens2math_extensions::Vector{<:Function}=Function[],
    ravens2math_passthrough::Dict{String,<:Vector{<:String}}=Dict{String,Vector{String}}(),
    make_pu_extensions::Vector{<:Function}=Function[],
    kwargs...
    )

    # @info "$(data)"

    @info "Converting CIM-RAVENS data model to MATHEMATICAL first to build JuMP model"

    data = transform_data_model_ravens(
            data;
            multinetwork=multinetwork,
            global_keys=global_keys,
            ravens2math_extensions=ravens2math_extensions,
            ravens2math_passthrough=ravens2math_passthrough,
            make_pu_extensions=make_pu_extensions,
        )

    @info "$(data["load"])"

    DEFINIDOEN_instantiate_mc_model

    return _IM.instantiate_model(
        data,
        model_type,
        build_method,
        ref_add_core!,
        union(_pmd_math_global_keys, global_keys),
        pmd_it_sym;
        ref_extensions=ref_extensions,
        kwargs...
    )
end

"""
    solve_mc_model(
        data::Dict{String,<:Any},
        model_type::Type,
        optimizer,
        build_mc::Function;
        ref_extensions::Vector{<:Function}=Function[],
        multinetwork::Bool=false,
        global_keys::Set{String}=Set{String}(),
        eng2math_extensions::Vector{<:Function}=Function[],
        eng2math_passthrough::Dict{String,<:Vector{<:String}}=Dict{String,Vector{String}}(),
        make_si::Bool=!get(data, "per_unit", false),
        make_si_extensions::Vector{<:Function}=Function[],
        dimensionalize_math_extensions::Dict{String,Dict{String,Vector{String}}}=Dict{String,Dict{String,Vector{String}}}(),
        kwargs...
    )::Dict{String,Any}

Takes data in either the ENGINEERING or MATHEMATICAL model, a model type (_e.g._, [`ACRUPowerModel`](@ref ACRUPowerModel)),
and model builder function (_e.g._, [`build_mc_opf`](@ref build_mc_opf)), and returns a solution in the original data model
defined by `data`.

If `make_si` is false, data will remain in per-unit.

For an explanation of `multinetwork` and `global_keys`, see [`make_multinetwork`](@ref make_multinetwork)

For an explanation of `eng2math_extensions` and `eng2math_passthrough`, see [`transform_data_model`](@ref transform_data_model)

For an explanation of `make_pu_extensions`, see [`make_per_unit!`](@ref make_per_unit!)

For an explanation of `ref_extensions`, see [`instantiate_mc_model`](@ref instantiate_mc_model)

For an explanation of `map_math2eng_extensions`, `make_si`, `make_si_extensions`, and `dimensionalize_math_extensions`, see [`solution_make_si`](@ref solution_make_si)
"""
function solve_mc_model(
    data::Dict{String,<:Any},
    model_type::Type,
    optimizer,
    build_mc::Function;
    ref_extensions::Vector{<:Function}=Function[],
    multinetwork::Bool=false,
    global_keys::Set{String}=Set{String}(),
    eng2math_extensions::Vector{<:Function}=Function[],
    eng2math_passthrough::Dict{String,<:Vector{<:String}}=Dict{String,Vector{String}}(),
    make_pu_extensions::Vector{<:Function}=Function[],
    map_math2eng_extensions::Dict{String,<:Function}=Dict{String,Function}(),
    make_si::Bool=!get(data, "per_unit", false),
    make_si_extensions::Vector{<:Function}=Function[],
    dimensionalize_math_extensions::Dict{String,Dict{String,Vector{String}}}=Dict{String,Dict{String,Vector{String}}}(),
    kwargs...
    )::Dict{String,Any}


    if iseng(data)
        data_math = transform_data_model(
            data;
            multinetwork=multinetwork,
            eng2math_extensions=eng2math_extensions,
            eng2math_passthrough=eng2math_passthrough,
            make_pu_extensions=make_pu_extensions,
            global_keys=global_keys,
        )

        result = _solve_mc_model(
            data_math,
            model_type,
            optimizer,
            build_mc;
            ref_extensions=ref_extensions,
            multinetwork=multinetwork,
            global_keys=global_keys,
            kwargs...
        )

        result["solution"] = transform_solution(
            result["solution"],
            data_math;
            map_math2eng_extensions=map_math2eng_extensions,
            make_si=make_si,
            make_si_extensions=make_si_extensions,
            dimensionalize_math_extensions=dimensionalize_math_extensions
        )
    elseif ismath(data)
        result = _solve_mc_model(
            data,
            model_type,
            optimizer,
            build_mc;
            ref_extensions=ref_extensions,
            multinetwork=multinetwork,
            eng2math_extensions=eng2math_extensions,
            eng2math_passthrough=eng2math_passthrough,
            global_keys=global_keys,
            kwargs...
        )
    else
        error("unrecognized data model format '$(get(data, "data_model", missing))'")
    end

    return result
end


"""
    solve_mc_model(
        file::String,
        model_type::Type,
        optimizer,
        build_mc::Function;
        dss2eng_extensions::Vector{<:Function}=Function[],
        multinetwork::Bool=false,
        global_keys::Set{String}=Set{String}(),
        kwargs...
    )::Dict{String,Any}

Given a `file::String`, data will be parsed automatically from the file.

See [`solve_mc_model`](@ref solve_mc_model) for detailed explanation of function arguments.
"""
function solve_mc_model(
    file::String,
    model_type::Type,
    optimizer,
    build_mc::Function;
    dss2eng_extensions::Vector{<:Function}=Function[],
    multinetwork::Bool=false,
    global_keys::Set{String}=Set{String}(),
    kwargs...
    )::Dict{String,Any}

    return solve_mc_model(
        parse_file(file; dss2eng_extensions=dss2eng_extensions, multinetwork=multinetwork, global_keys=global_keys),
        model_type,
        optimizer,
        build_mc;
        multinetwork=multinetwork,
        global_keys=global_keys,
        kwargs...
    )
end
