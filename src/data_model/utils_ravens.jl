
"initializes the base math object of any type"
function _init_math_obj_ravens(obj_type::String, eng_id::Any, eng_obj::Dict{String,<:Any}, index::Int; pass_props::Vector{String}=String[])::Dict{String,Any}
    math_obj = Dict{String,Any}(
        "name" => "$eng_id",
        "source_id" => "$obj_type.$eng_id"
    )

    math_obj["index"] = index

    return math_obj
end


"converts impendance in Ohm/m by multiplying by length"
function _impedance_conversion_ravens(data_eng::Dict{String,Any}, eng_obj::Dict{String,Any}, key::String)

    _conductor_count =  data_eng["PerLengthPhaseImpedance.conductorCount"]
    _impedance_matrix = zeros(Float64, _conductor_count, _conductor_count)

    _index = 1
    for i in 1:_conductor_count
        for j in 1:i
            _impedance_matrix[i, j] =  get(data_eng["PerLengthPhaseImpedance.PhaseImpedanceData"][_index], key, 0.0)
            _impedance_matrix[j, i] = get(data_eng["PerLengthPhaseImpedance.PhaseImpedanceData"][_index], key, 0.0)
            _index += 1
        end
    end

    return _impedance_matrix .* get(eng_obj, "Conductor.length", 1.0)
end


"converts admittance by multiplying by 2πωl"
function _admittance_conversion_ravens(data_eng::Dict{String,<:Any}, eng_obj::Dict{String,<:Any}, key::String; freq::Float64=60.0)

    _conductor_count =  data_eng["PerLengthPhaseImpedance.conductorCount"]
    _admittance_matrix = zeros(Float64, _conductor_count, _conductor_count)

    _index = 1
    for i in 1:_conductor_count
        for j in 1:i
            _admittance_matrix[i, j] =  get(data_eng["PerLengthPhaseImpedance.PhaseImpedanceData"][_index], key, 0.0)
            _admittance_matrix[j, i] = get(data_eng["PerLengthPhaseImpedance.PhaseImpedanceData"][_index], key, 0.0)
            _index += 1
        end
    end

    return _admittance_matrix .* get(eng_obj, "Conductor.length", 1.0) .* freq ./ 1e2 # divide by 2 to get both sides _to and _fr
end
