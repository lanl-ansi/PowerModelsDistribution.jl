"""
    parse_file(io)

Parses the IOStream of a file into a PowerModelsDistribution data structure.
"""
function parse_file(io::IO, filetype::AbstractString="json"; data_model::DataModel=ENGINEERING, import_all::Bool=false, bank_transformers::Bool=true, transformations::Vector{<:Function}=Vector{Function}([]))::Dict{String,Any}
    if filetype == "dss"
        data_eng = PowerModelsDistribution.parse_opendss(io; import_all=import_all, bank_transformers=bank_transformers)

        for transformation in transformations
            transformation(data_eng)
        end

        if data_model == MATHEMATICAL
            return transform_data_model(data_eng; make_pu=true)
        else
            return data_eng
        end
    elseif filetype == "json"
        pmd_data = parse_json(io; validate=false)

        if pmd_data["data_model"] != data_model && data_model == ENGINEERING
            return transform_data_model(pmd_data)
        else
            return pmd_data
        end
    else
        Memento.error(_LOGGER, "only .dss and .json files are supported")
    end
end


""
function parse_file(file::String; kwargs...)::Dict{String,Any}
    data = open(file) do io
        parse_file(io, split(lowercase(file), '.')[end]; kwargs...)
    end

    return data
end


"transforms model between engineering (high-level) and mathematical (low-level) models"
function transform_data_model(data::Dict{String,<:Any}; kron_reduced::Bool=true, make_pu::Bool=true)::Dict{String,Any}
    current_data_model = get(data, "data_model", MATHEMATICAL)

    if current_data_model == ENGINEERING
        data_math = _map_eng2math(data; kron_reduced=kron_reduced)

        correct_network_data!(data_math; make_pu=make_pu)

        return data_math
    elseif current_data_model == MATHEMATICAL
        Memento.warn(_LOGGER, "A MATHEMATICAL data model cannot be converted back to an ENGINEERING data model, irreversible transformations have been made")
        return data
    else
        Memento.warn(_LOGGER, "Data model '$current_data_model' is not recognized, no model type transformation performed")
        return data
    end
end


""
function correct_network_data!(data::Dict{String,Any}; make_pu::Bool=true)
    if get(data, "data_model", MATHEMATICAL) == ENGINEERING
        check_eng_data_model(data)
    else
        if make_pu
            make_per_unit!(data)

            _PM.check_connectivity(data)
            _PM.correct_transformer_parameters!(data)
            _PM.correct_voltage_angle_differences!(data)
            _PM.correct_thermal_limits!(data)
            _PM.correct_branch_directions!(data)
            _PM.check_branch_loops(data)
            _PM.correct_bus_types!(data)
            _PM.correct_dcline_limits!(data)
            _PM.correct_cost_functions!(data)
            _PM.standardize_cost_terms!(data)
        end
    end
end
