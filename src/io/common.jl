"""
    parse_file(io)

Parses the IOStream of a file into a Three-Phase PowerModels data structure.
"""
function parse_file(io::IO, filetype::AbstractString="json"; data_model::String="mathematical", import_all::Bool=false, bank_transformers::Bool=true, lossless::Bool=false, use_dss_bounds::Bool=true)::Dict{String,Any}
    if filetype == "dss"
        data_eng = PowerModelsDistribution.parse_opendss(io; import_all=import_all, bank_transformers=bank_transformers)

        if data_model == "mathematical"
            return transform_data_model(data_eng; make_pu=true, lossless=lossless, use_dss_bounds=use_dss_bounds)
        else
            return data_eng
        end
    elseif filetype == "json"
        pmd_data = parse_json(io; validate=false)

        if get(pmd_data, "data_model", "mathematical") != data_model
            return transform_data_model(pmd_data; lossless=lossless, use_dss_bounds=use_dss_bounds)
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
function transform_data_model(data::Dict{String,<:Any}; kron_reduced::Bool=true, make_pu::Bool=false, lossless::Bool=false, use_dss_bounds::Bool=true)::Dict{String,Any}
    current_data_model = get(data, "data_model", "mathematical")

    if current_data_model == "engineering"
        data_math = _map_eng2math(data; kron_reduced=kron_reduced, lossless=lossless, use_dss_bounds=use_dss_bounds)

        correct_network_data!(data_math; make_pu=make_pu)

        return data_math
    elseif current_data_model == "mathematical"
        data_eng = _map_math2eng(data)

        correct_network_data!(data_eng; make_pu=make_pu)
    else
        @warn "Data model '$current_data_model' is not recognized, no model type transformation performed"
        return data
    end
end


""
function correct_network_data!(data::Dict{String,Any}; make_pu::Bool=true)
    if get(data, "data_model", "mathematical") == "engineering"
        check_eng_data_model(data)
        if make_pu
            make_per_unit!(data)
        end
    else
        if make_pu
            make_per_unit!(data)

            _PMs.check_connectivity(data)
            _PMs.correct_transformer_parameters!(data)
            _PMs.correct_voltage_angle_differences!(data)
            _PMs.correct_thermal_limits!(data)
            _PMs.correct_branch_directions!(data)
            _PMs.check_branch_loops(data)
            _PMs.correct_bus_types!(data)
            _PMs.correct_dcline_limits!(data)
            _PMs.correct_cost_functions!(data)
            _PMs.standardize_cost_terms!(data)
        end
    end
end
