const _phasecode_map = Dict(
    "PhaseCode.ABC" => [1, 2, 3],
    "PhaseCode.AB" => [1, 2],
    "PhaseCode.AC" => [1, 3],
    "PhaseCode.BC" => [2, 3],
    "PhaseCode.A" => [1],
    "PhaseCode.B" => [2],
    "PhaseCode.C" => [3]
)

_phase_map = Dict(
    "SinglePhaseKind.A" => 1,
    "SinglePhaseKind.B" => 2,
    "SinglePhaseKind.C" => 3,
    "SinglePhaseKind.N" => 4
)

const _multipliers_map = Dict(
    "m" => 1e-3,
    "c" => 1e-2,
    "d" => 1e-1,
    "da" => 1e1,
    "h" => 1e2,
    "k" => 1e3,
    "M" => 1e6,
    "G" => 1e9,
    "T" => 1e12,
    "P" => 1e15,
    "E" => 1e18,
    "Z" => 1e21
)


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
function _impedance_conversion_ravens(eng_obj::Dict{String,<:Any}, vals::Matrix{Float64})
    return vals .* get(eng_obj, "Conductor.length", 1.0)
end


"converts admittance by multiplying by 2πωl"
function _admittance_conversion_ravens(eng_obj::Dict{String,<:Any}, vals::Matrix{Float64})
    2.0 .* pi .* vals .* get(eng_obj, "Conductor.length", 1.0) ./ 1e9
end


"converts impendance in Ohm/m by multiplying by length"
function _impedance_conversion_ravens(data_eng::Dict{String,Any}, eng_obj::Dict{String,Any}, key::String)

    _conductor_count =  data_eng["PerLengthPhaseImpedance.conductorCount"]
    _impedance_matrix = zeros(Float64, _conductor_count, _conductor_count)

    for obj in data_eng["PerLengthPhaseImpedance.PhaseImpedanceData"]
        row = obj["PhaseImpedanceData.row"]
        col = obj["PhaseImpedanceData.column"]
        value = get(obj, key, 0.0)
        _impedance_matrix[row, col] = value
        _impedance_matrix[col, row] = value
    end

    return _impedance_matrix .* get(eng_obj, "Conductor.length", 1.0)
end


"converts impendance in Ohm/m by multiplying by length"
function _impedance_conversion_ravens_energy_source(data_eng::Dict{String,Any}, eng_obj::Dict{String,Any}, key1::String, key2::String)
    # Default energy sources considered 3 phases
    nphases = 3
    _impedance_matrix = zeros(Float64, nphases, nphases)

    z = get(eng_obj, key1, 0.0)
    z0 = get(eng_obj, key2, 0.0)

    for i in 1:nphases
        for j in 1:i
            if(i==j)
                _impedance_matrix[i, j] =  z + ((z0 - z)/3)
            else
                _impedance_matrix[i, j] = (z0 - z)/3
                _impedance_matrix[j, i] = (z0 - z)/3
            end
        end
    end

    return _impedance_matrix .* get(eng_obj, "Conductor.length", 1.0)
end


"converts admittance by multiplying by 2πωl"
function _admittance_conversion_ravens(data_eng::Dict{String,<:Any}, eng_obj::Dict{String,<:Any}, key::String)

    _conductor_count =  data_eng["PerLengthPhaseImpedance.conductorCount"]
    _admittance_matrix = zeros(Float64, _conductor_count, _conductor_count)

    for obj in data_eng["PerLengthPhaseImpedance.PhaseImpedanceData"]
        row = obj["PhaseImpedanceData.row"]
        col = obj["PhaseImpedanceData.column"]
        value = get(obj, key, 0.0)
        _admittance_matrix[row, col] = value
        _admittance_matrix[col, row] = value
    end

    return _admittance_matrix .* get(eng_obj, "Conductor.length", 1.0) ./ 2 # divide by 2 to get both sides _to and _fr
end

"extracts the name from a ravens reference string"
function _extract_name(element)

    name = replace(split(element, "::")[2], "'" => "")
    return name
end


"calculates the shunt admittance matrix based on terminals and b or g"
function _calc_shunt_admittance_matrix(terminals, b)

    N = length(terminals)
    _shunt_matrix = b* Matrix(LinearAlgebra.I, N, N)
    return _shunt_matrix

end
