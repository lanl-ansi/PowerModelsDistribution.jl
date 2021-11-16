"""
    parse_file(
        io::IO,
        filetype::Union{AbstractString,Missing}=missing;
        data_model::DataModel=ENGINEERING,
        import_all::Bool=false,
        bank_transformers::Bool=true,
        transformations::Vector{<:Any}=[],
        dss2eng_extensions::Vector{<:Function}=Function[],
        eng2math_extensions::Vector{<:Function}=Function[],
        eng2math_passthrough::Dict{String,Vector{String}}=Dict{String,Vector{String}}(),
        make_pu_extensions::Vector{<:Function}=Function[],
        make_pu::Bool=true,
        multinetwork::Bool=false,
        global_keys::Set{String}=Set{String}(),
        kron_reduce::Bool=true,
        phase_project::Bool=false,
        time_series::String="daily"
    )::Dict{String,Any}

Parses the IOStream of a file into a PowerModelsDistribution data structure

If `filetype` is missing, `parse_file` will attempt to detect the filetype, but this may fail, and it is
advised to pass the filetype if it is known.

If `data_model` is `MATHEMATICAL`, the data model type will be automatically transformed via
[`transform_data_model`](@ref transform_data_model).

For explanation of `import_all`, `bank_transformers`, and `time_series`, see [`parse_opendss`](@ref parse_opendss)

For explanation of `dss2eng_extensions`, see [`parse_opendss`](@ref parse_opendss)

For explanation of `kron_reduce`, see [`apply_kron_reduction!`](@ref apply_kron_reduction!)

For explanation of `phase_project`, see [`apply_phase_projection!`](@ref apply_phase_projection!)

For explanation of `multinetwork` and `global_keys`, see [`make_multinetwork`](@ref make_multinetwork) and [`transform_data_model`](@ref transform_data_model)

For explanation of `eng2math_extensions` and `eng2math_passthrough`, see [`transform_data_model`](@ref transform_data_model)

For explanation of `make_pu` and `make_pu_extensions`, see [`make_per_unit!`](@ref make_per_unit!).
"""
function parse_file(
    io::IO,
    filetype::Union{AbstractString,Missing}=missing;
    data_model::DataModel=ENGINEERING,
    import_all::Bool=false,
    bank_transformers::Bool=true,
    transformations::Vector{<:Any}=[],
    dss2eng_extensions::Vector{<:Function}=Function[],
    eng2math_extensions::Vector{<:Function}=Function[],
    eng2math_passthrough::Dict{String,Vector{String}}=Dict{String,Vector{String}}(),
    make_pu_extensions::Vector{<:Function}=Function[],
    make_pu::Bool=true,
    multinetwork::Bool=false,
    global_keys::Set{String}=Set{String}(),
    kron_reduce::Bool=true,
    phase_project::Bool=false,
    time_series::String="daily"
    )::Dict{String,Any}

    pmd_data = Dict{String,Any}()
    if ismissing(filetype) || filetype == "json"
        try
            pmd_data = parse_json(io)
            filetype = "json"
        catch err
            filetype = "dss"
        end
    end

    if filetype == "dss"
        pmd_data = parse_opendss(io;
            import_all=import_all,
            bank_transformers=bank_transformers,
            time_series=time_series,
            dss2eng_extensions=dss2eng_extensions,
        )

        for transform! in transformations
            @assert isa(transform!, Function) || isa(transform!, Tuple{<:Function,Vararg{Pair{String,<:Any}}})

            if isa(transform!, Tuple)
                transform![1](pmd_data; [Symbol(k)=>v for (k,v) in transform![2:end]]...)
            else
                transform!(pmd_data)
            end
        end

        if multinetwork
            pmd_data = make_multinetwork(pmd_data; global_keys=global_keys)
        end

        if data_model == MATHEMATICAL
            pmd_data = transform_data_model(
                pmd_data;
                make_pu=make_pu,
                make_pu_extensions=make_pu_extensions,
                kron_reduce=kron_reduce,
                phase_project=phase_project,
                multinetwork=multinetwork,
                global_keys=global_keys,
                eng2math_extensions=eng2math_extensions,
                eng2math_passthrough=eng2math_passthrough,
            )
        end
    elseif filetype == "json"
        if multinetwork && !ismultinetwork(pmd_data)
            pmd_data = make_multinetwork(pmd_data; global_keys=global_keys)
        end

        if data_model == MATHEMATICAL && !ismath(pmd_data)
            pmd_data = transform_data_model(pmd_data;
                make_pu=make_pu,
                make_pu_extensions=make_pu_extensions,
                kron_reduce=kron_reduce,
                phase_project=phase_project,
                multinetwork=multinetwork,
                global_keys=global_keys,
                eng2math_extensions=eng2math_extensions,
                eng2math_passthrough=eng2math_passthrough,
            )
        end
    else
        error("only .dss and .json files are supported")
    end

    return pmd_data
end


"""
    parse_file(file::String; kwargs...)::Dict{String,Any}

Loads file into IOStream and passes it onto [`parse_file`](@ref parse_file)
"""
function parse_file(file::String; kwargs...)::Dict{String,Any}
    data = open(file) do io
        parse_file(io, split(lowercase(file), '.')[end]; kwargs...)
    end

    return data
end


"""
    correct_network_data!(data::Dict{String,Any}; make_pu::Bool=true, make_pu_extensions::Vector{<:Function}=Function[])

Makes corrections and performs checks on network data structure in either ENGINEERING or MATHEMATICAL formats, and converts
to per-unit if data a is MATHEMATICAL data model.

If `make_pu` is false, converting to per-unit will be skipped.

# Custom per-unit transformations

See [`make_per_unit!`](@ref make_per_unit!)
"""
function correct_network_data!(data::Dict{String,Any}; make_pu::Bool=true, make_pu_extensions::Vector{<:Function}=Function[])
    if iseng(data)
        check_eng_data_model(data)
    elseif ismath(data)
        check_connectivity(data)

        correct_branch_directions!(data)
        check_branch_loops(data)

        correct_bus_types!(data)

        propagate_network_topology!(data)

        if make_pu
            make_per_unit!(data; make_pu_extensions=make_pu_extensions)

            correct_mc_voltage_angle_differences!(data)
            correct_mc_thermal_limits!(data)

            correct_cost_functions!(data)
            standardize_cost_terms!(data)
        end
    end
end
