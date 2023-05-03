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
    data_model::Union{DataModel,Type{T}}=EngineeringModel{NetworkModel},
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
    ) where T <: DistributionModel

    if isa(data_model, DataModel)
        if data_model == DSS
            data_model = DssModel
        elseif data_model == ENGINEERING
            if multinetwork
                data_model = EngineeringModel{MultinetworkModel}
            else
                data_model = EngineeringModel{NetworkModel}
            end
        elseif data_model == MATHEMATICAL
            if multinetwork
                data_model = MathematicalModel{MultinetworkModel}
            else
                data_model = MathematicalModel{NetworkModel}
            end
        else
            throw(ErrorException("Unsupported data model $(data_model)"))
        end
    end

    data = missing
    if ismissing(filetype) || filetype == "json"
        try
            data = parse_json(io)
            filetype = "json"
        catch err
            if isa(err, ArgumentError)
                filetype = "dss"
            else
                throw(err)
            end
        end
    end

    if filetype == "dss"
        data = parse_raw_dss(io)
    elseif filetype != "json"
        error("only .dss and .json files are supported")
    end

    if !ismissing(data)
        transform_data_model(
            data_model,
            data;
            dss2eng_extensions=dss2eng_extensions,
            kron_reduce=kron_reduce,
            phase_project=phase_project,
            time_series=time_series,
            import_all=import_all,
            transformations=transformations,
            bank_transformers=bank_transformers,
            eng2math_extensions=eng2math_extensions,
            eng2math_passthrough=eng2math_passthrough,
            make_pu=make_pu,
            make_pu_extensions=make_pu_extensions,
            global_keys=global_keys,
        )
    else
        error("Filetype is not recognized, no data loaded")
    end
end


function read(io::IO, ::Type{DssModel}, ::Type{EngineeringModel{NetworkModel}})::DistributionModel
end


"""
    parse_file(file::String; kwargs...)::Dict{String,Any}

Loads file into IOStream and passes it onto [`parse_file`](@ref parse_file)
"""
function parse_file(file::String; kwargs...)
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

correct_network_data!(eng::EngineeringDataModel) = nothing
correct_network_data!(math::MathematicalModel) = nothing
