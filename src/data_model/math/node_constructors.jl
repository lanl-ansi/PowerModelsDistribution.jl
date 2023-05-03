function create_math_object(::Type{T}, eng_obj::EngBus, id::Int, bus_lookup::Dict{String,Int}, eng::EngineeringModel{NetworkModel})::Tuple{T,Vector{MathShunt}} where T <: MathBus
    grounded_perfect, shunts = _convert_grounding(eng_obj.terminals, eng_obj.grounded, eng_obj.rg, eng_obj.xg)

    vmin, vmax = _get_tight_absolute_voltage_magnitude_bounds(eng_obj)
    vm_pair_lb, vm_pair_ub  = _get_tight_pairwise_voltage_magnitude_bounds(eng_obj)

    to_sh = MathShunt[]
    for (i,(sh_connections, sh_y)) in enumerate(shunts)
        push!(to_sh, EngShuntObj(;
            index = 1,
            shunt_bus = id,
            connections = sh_connections,
            gs = real.(sh_y),
            bs = imag.(sh_y),
            source_id = "$(typeof(T)).$(eng_obj.name).ground.$i"
        ))
    end

    math_obj = T(;
        index = id,
        bus_i = id,
        terminals = eng_obj.terminals,
        grounded = grounded_perfect,
        vm = eng_obj.vm,
        va = eng_obj.va,
        vmin = vmin,
        vmax = vmax,
        vm_pair_lb = vm_pair_lb,
        vm_pair_ub = vm_pair_ub,
        bus_type = eng_obj.status == DISABLED ? 4 : 1,
        source_id = "$(typeof(eng_obj)).$(eng_obj.name)",
    )

    return math_obj, to_sh
end

function create_math_object(::Type{T}, eng_obj::EngLoad, id::Int, bus_lookup::Dict{String,Int}, eng::EngineeringModel{NetworkModel})::Union{Vector{T},T} where T <: MathLoad
    if eng_obj.model == ZIP
        math_objs = T[]
        for (idx,l) in enumerate([IMPEDANCE,CURRENT,POWER])
            push!(math_objs, T(;
                model = l,
                connections = eng_obj.connections,
                configuration = eng_obj.configuration,
                status = Int(eng_obj.status),
                pd = eng_obj.pd_nom*eng_obj.zipv[idx],
                qd = eng_obj.qd_nom*eng_obj.zipv[idx],
                load_bus = bus_lookup[eng_obj.bus],
                dispatchable = Int(eng_obj.dispatchable),
                source_id = "$(typeof(eng_obj)).$(eng_obj.name)",
                index = id+(idx-1),
                vnom_kv = eng_obj.vm_nom,
            ))
        end

        return math_objs
    else
        T(;
            index = id,
            load_bus = bus_lookup[eng_obj.bus],
            connections = eng_obj.connections,
            configuration = eng_obj.configuration,
            model = eng_obj.model,
            pd = eng_obj.pd_nom,
            qd = eng_obj.qd_nom,
            vnom_kv = eng_obj.vm_nom,
            dispatchable = Int(eng_obj.dispatchable),
            status = Int(eng_obj.status),
            source_id = "$(typeof(eng_obj)).$(eng_obj.name)",
            dss = eng_obj.dss
        )
    end
end

function create_math_object(::Type{T}, eng_obj::EngShunt, id::Int, bus_lookup::Dict{String,Int}, eng::EngineeringModel{NetworkModel})::T where T <: MathShunt
    math_obj = T(;
        index = id,
        shunt_bus = bus_lookup[eng_obj.bus],
        connections = eng_obj.connections,
        source_id = "$(typeof(T)).$(eng_obj.name)",
        gs = eng_obj.gs,
        bs = eng_obj.bs,
        # controls = create_math_object(MathShuntControlsObj, eng_obj.controls, 0, bus_lookup), # TODO
        status = Int(eng_obj.status),
        dispatchable = Int(eng_obj.dispatchable),
        dss = eng_obj.dss
    )
end

function create_math_object(::Type{T}, eng_obj::EngGenerator, id::Int, bus_lookup::Dict{String,Int}, eng::EngineeringModel{NetworkModel})::T where T <: MathGenerator
    math_obj = T(;
        gen_bus = bus_lookup[eng_obj.bus],
        index = id,
        gen_status = Int(eng_obj.status),
        control_mode = Int(eng_obj.control_mode),
        pmin = eng_obj.pg_lb,
        pmax = eng_obj.pg_ub,
        qmin = eng_obj.qg_lb,
        qmax = eng_obj.qg_ub,
        pg = eng_obj.pg,
        qg = eng_obj.qg,
        vg = eng_obj.vg,
        configuration = eng_obj.configuration,
        model = eng_obj.cost_pg_model,
        startup = 0.0,
        shutdown = 0.0,
        cost = eng_obj.cost_pg_parameters,
        ncost = length(eng_obj.cost_pg_parameters),
        source_id = "$(typeof(T)).$(eng_obj.name)",
        dss = eng_obj.dss
    )
end

function create_math_object(::Type{T}, eng_obj::EngSolar, id::Int, bus_lookup::Dict{String,Int}, eng::EngineeringModel{NetworkModel})::T where T <: MathGenerator
    math_obj = T(;
        gen_bus = bus_lookup[eng_obj.bus],
        index = id,
        gen_status = Int(eng_obj.status),
        control_mode = Int(eng_obj.control_mode),
        pmin = eng_obj.pg_lb,
        pmax = eng_obj.pg_ub,
        qmin = eng_obj.qg_lb,
        qmax = eng_obj.qg_ub,
        pg = eng_obj.pg,
        qg = eng_obj.qg,
        vg = eng_obj.vg,
        configuration = eng_obj.configuration,
        model = eng_obj.cost_pg_model,
        startup = 0.0,
        shutdown = 0.0,
        cost = eng_obj.cost_pg_parameters,
        ncost = length(eng_obj.cost_pg_parameters),
        source_id = "$(typeof(T)).$(eng_obj.name)",
        dss = eng_obj.dss
    )
end

function create_math_object(::Type{T}, eng_obj::EngVoltageSource, id::Int, bus_lookup::Dict{String,Int}, eng::EngineeringModel{NetworkModel})::Union{T,Tuple{T,MathBus,MathBranch}} where T <: MathGenerator
    math_obj = T(;
        index = id,
        gen_bus = bus_lookup[eng_obj.bus],
        gen_status = Int(eng_obj.status),
        control_mode = Int(eng_obj.control_mode),
        pmin = eng_obj.pg_lb,
        pmax = eng_obj.pg_ub,
        qmin = eng_obj.qg_lb,
        qmax = eng_obj.qg_ub,
        pg = eng_obj.pg,
        qg = eng_obj.qg,
        vg = eng_obj.vg,
        configuration = eng_obj.configuration,
        model = eng_obj.cost_pg_model,
        startup = 0.0,
        shutdown = 0.0,
        cost = eng_obj.cost_pg_parameters,
        ncost = length(eng_obj.cost_pg_parameters),
        source_id = "$(typeof(T)).$(eng_obj.name)",
        dss = eng_obj.dss
    )

    if !all(isapprox.(get(eng_obj, "rs", zeros(1, 1)), 0)) && !all(isapprox.(get(eng_obj, "xs", zeros(1, 1)), 0))
        bus_obj = MathBusObj(;
            index = length(bus_lookup)+1,
            bus_i = length(bus_lookup)+1,
            terminals = eng_obj.connections,
            grounded = fill(false, length(eng_obj.connections)),
            bus_type = Int(eng_obj.status) == 0 ? 4 : math_obj.control_mode == Int(ISOCHRONOUS) ? 3 : 2,
            vm = eng_obj.vm,
            va = eng_obj.va,
            vmin = eng_obj.vm_lb,
            vmax = eng_obj.vm_ub,
            source_id = "$(typeof(T)).$(eng_obj.name)",
        )

        branch_obj = MathBranchObj(;
            index = 0,
            f_bus = math_obj.gen_bus,
            t_bus = bus_obj.bus_i,
            f_connections = eng_obj.connections,
            t_connections = eng_obj.connections,
            angmin = fill(-10.0, length(eng_obj.connections)),
            angmax = fill( 10.0, length(eng_obj.connections)),
            br_status = math_obj.gen_status,
            br_r = eng_obj.rs,
            br_x = eng_obj.xs,
        )

        math_obj.gen_bus = bus_obj.bus_i

        return (math_obj, bus_obj, branch_obj)
    end

    return math_obj
end

function create_math_object(::Type{T}, eng_obj::EngStorage, id::Int, bus_lookup::Dict{String,Int}, eng::EngineeringModel{NetworkModel})::T where T <: MathStorage
    math_obj = T(;
        storage_bus = bus_lookup[eng_obj.bus],
        index = id,
        configuration = eng_obj.configuration,
        energy = eng_obj.energy,
        energy_rating = eng_obj.energy_ub,
        charge_rating = eng_obj.charge_ub,
        discharge_rating = eng_obj.discharge_ub,
        charge_efficiency = eng_obj.charge_efficiency / 100.0,
        discharge_efficiency = eng_obj.discharge_efficiency / 100.0,
        thermal_rating = eng_obj.sm_ub,
        qmin = eng_obj.qs_lb,
        qmax = eng_obj.qs_ub,
        r = eng_obj.rs,
        x = eng_obj.xs,
        p_loss = eng_obj.pex,
        q_loss = eng_obj.qex,
        ps = eng_obj.ps,
        qs = eng_obj.qs,
        control_mode = Int(eng_obj.control_mode),
        status = Int(eng_obj.status),
        dispatchable = Int(eng_obj.dispatchable),
        source_id = "$(typeof(T)).$(eng_obj.name)",
        dss = eng_obj.dss
    )
end
