function create_math_object(::Type{T}, eng_obj::EngLine, id::Int, bus_lookup::Dict{String,Int}; is_per_unit::Bool=false) where T <: MathBranch
end

function create_math_object(::Type{T}, eng_obj::EngSwitch, id::Int, bus_lookup::Dict{String,Int}; is_per_unit::Bool=false) where T <: MathSwitch
end

function create_math_object(::Type{T}, eng_obj::EngTransformer, id::Int, bus_lookup::Dict{String,Int}; is_per_unit::Bool=false) where T <: MathTransformer
end
