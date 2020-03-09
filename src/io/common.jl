"""
    parse_file(io)

Parses the IOStream of a file into a Three-Phase PowerModels data structure.
"""
function parse_file(io::IO; data_model::String="mathematical", import_all::Bool=false, filetype::AbstractString="json", bank_transformers::Bool=true)
    if filetype == "m"
        pmd_data = PowerModelsDistribution.parse_matlab(io)
    elseif filetype == "dss"
        Memento.warn(_LOGGER, "Not all OpenDSS features are supported, currently only minimal support for lines, loads, generators, and capacitors as shunts. Transformers and reactors as transformer branches are included, but value translation is not fully supported.")
        pmd_data = PowerModelsDistribution.parse_opendss(io; import_all=import_all, bank_transformers=bank_transformers)
    elseif filetype == "json"
        pmd_data = PowerModels.parse_json(io; validate=false)
    else
        Memento.error(_LOGGER, "only .m and .dss files are supported")
    end

    if data_model == "mathematical"
        pmd_data = transform_data_model(pmd_data)

        correct_network_data!(pmd_data)
    end

    return pmd_data
end


""
function parse_file(file::String; kwargs...)
    pmd_data = open(file) do io
        parse_file(io; filetype=split(lowercase(file), '.')[end], kwargs...)
    end

    return pmd_data
end


"transforms model between engineering (high-level) and mathematical (low-level) models"
function transform_data_model(data::Dict{<:Any,<:Any}; kron_reduced::Bool=true)
    current_data_model = get(data, "data_model", "mathematical")

    if current_data_model == "engineering"
        out = _map_eng2math(data, kron_reduced=kron_reduced)

        bus_indexed_id = string(out["bus_lookup"][data["settings"]["set_vbase_bus"]])
        vbases = Dict(bus_indexed_id=>data["settings"]["set_vbase_val"])
        sbase = data["settings"]["set_sbase_val"]

        make_per_unit!(out, vbases=vbases, sbase=sbase, v_var_scalar=data["settings"]["v_var_scalar"])
        return out
    elseif current_data_model == "mathematical"
        return _map_math2eng!(data)
    else
        @warn "Data model '$current_data_model' is not recognized, no transformation performed"
        return data
    end
end


""
function correct_network_data!(data::Dict{String,Any})
    #_PMs.make_per_unit!(data)

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
