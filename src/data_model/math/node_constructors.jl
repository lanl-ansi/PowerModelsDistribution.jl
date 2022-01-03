function create_math_object(::Type{T}, eng_obj::EngBus, id::Int) where T <: MathBus
    grounded_perfect, shunts = _convert_grounding(eng_obj.terminals, eng_obj.grounded, eng_obj.rg, eng_obj.xg)

    T(;
        id = id,
        terminals = eng_obj.terminals,
        grounded = grounded_perfect,
        vm = eng_obj.vm,
        va = eng_obj.va,
        vm_lb = eng_obj.vm_lb,
        vm_ub = eng_obj.vm_ub,
        vm_pair_lb = eng_obj.vm_pair_lb,
        vm_pair_ub = eng_obj.vm_pair_ub,
        status = eng_obj.status,
        source_id = "EngBus.$(eng_obj.name)",
    )
end

function create_math_object(::Type{T}, eng_obj::EngLoad, id::Int, bus_lookup::Dict{String,Int}) where T <: MathLoad
    T(;
        id = id,
        bus = bus_lookup[eng_obj.name],
        connections = eng_obj.connections,
        configuration = eng_obj.configuration,
        model = eng_obj.model,
        pd = eng_obj.pd_nom,
        qd = eng_obj.qd_nom,
        vm_nom = eng_obj.vm_nom,
        dispatchable = eng_obj.dispatchable,
        status = eng_obj.status,
        source_id = "EngLoad.$(eng_obj.name)",
    )
end

function create_math_object(::Type{T}, eng_obj::EngShunt, id::Int, bus_lookup::Dict{String,Int}) where T <: MathShunt
end

function create_math_object(::Type{T}, eng_obj::EngGenerator, id::Int, bus_lookup::Dict{String,Int}) where T <: MathGenerator
end

function create_math_object(::Type{T}, eng_obj::EngSolar, id::Int, bus_lookup::Dict{String,Int}) where T <: MathGenerator
end

function create_math_object(::Type{T}, eng_obj::EngVoltageSource, id::Int, bus_lookup::Dict{String,Int}) where T <: MathGenerator
end

function create_math_object(::Type{T}, eng_obj::EngStorage, id::Int, bus_lookup::Dict{String,Int}) where T <: MathStorage
end
