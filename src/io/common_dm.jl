"""
    parse_file(io)

Parses the IOStream of a file into a Three-Phase PowerModels data structure.
"""
function parse_file_dm(io::IO; import_all::Bool=false, filetype::AbstractString="json", bank_transformers::Bool=true)
    if filetype == "m"
        pmd_data = PowerModelsDistribution.parse_matlab(io)
    elseif filetype == "dss"
        Memento.warn(_LOGGER, "Not all OpenDSS features are supported, currently only minimal support for lines, loads, generators, and capacitors as shunts. Transformers and reactors as transformer branches are included, but value translation is not fully supported.")
        pmd_data = PowerModelsDistribution.parse_opendss_dm(io; import_all=import_all, bank_transformers=bank_transformers)
    elseif filetype == "json"
        pmd_data = PowerModels.parse_json(io; validate=false)
    else
        Memento.error(_LOGGER, "only .m and .dss files are supported")
    end

    #correct_network_data!(pmd_data)

    return pmd_data
end


""
function parse_file_dm(file::String; data_model::String="engineering", kwargs...)
    pmd_data = open(file) do io
        parse_file_dm(io; filetype=split(lowercase(file), '.')[end], kwargs...)
    end

    if data_model == "mathematical"
        # to get the old indexing
        pmd_data_old = open(file) do io
            parse_file(io; filetype=split(lowercase(file), '.')[end], kwargs...)
        end
        index_presets = Dict(comp_type=>Dict(comp["name"]=>comp["index"] for (id, comp) in pmd_data_old[comp_type]) for comp_type in ["bus"])

        transform_data_model!(pmd_data; index_presets=index_presets)
    end

    return pmd_data
end


"transforms model between engineering (high-level) and mathematical (low-level) models"
function transform_data_model!(data::Dict{String,<:Any}; index_presets::Dict{String,<:Any}=Dict{String,Any}())
    if get(data, "data_model", "mathematical") == "engineering"
        data_model_map!(data)
        data_model_make_pu!(data)
        data_model_index!(data, index_presets=index_presets)
        data_model_make_compatible_v8!(data)
    else
        if haskey(data, "")
        else
            Memento.warn(_LOGGER, "Cannot transform mathematical model to engineering model, no mapping information available")
        end
    end
end
