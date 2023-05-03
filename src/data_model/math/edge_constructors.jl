function create_math_object(::Type{T}, eng_obj::EngLine, id::Int, bus_lookup::Dict{String,Int}, eng::EngineeringModel{NetworkModel}; is_per_unit::Bool=false)::T where T <: MathBranch
    math_obj = T(;
        index = id,
        f_bus = bus_lookup[eng_obj.f_bus],
        t_bus = bus_lookup[eng_obj.t_bus],
        f_connections = eng_obj.f_connections,
        t_connections = eng_obj.t_connections,
        br_r = eng_obj.rs .* eng_obj.length,
        br_x = eng_obj.xs .* eng_obj.length,
        g_fr = 2π .* eng.settings.base_frequency .* eng_obj.length .* eng_obj.g_fr ./ 1e9,
        g_to = 2π .* eng.settings.base_frequency .* eng_obj.length .* eng_obj.g_to ./ 1e9,
        b_fr = 2π .* eng.settings.base_frequency .* eng_obj.length .* eng_obj.b_fr ./ 1e9,
        b_to = 2π .* eng.settings.base_frequency .* eng_obj.length .* eng_obj.b_to ./ 1e9,
        angmin = eng_obj.vad_lb,
        angmax = eng_obj.vad_ub,
        c_rating_a = eng_obj.cm_ub,
        rate_a = eng_obj.sm_ub,
        is_per_unit=is_per_unit,
        br_status = Int(eng_obj.status),
        dss = eng_obj.dss
    )
end

function create_math_object(::Type{T}, eng_obj::EngSwitch, id::Int, bus_lookup::Dict{String,Int}, eng::EngineeringModel{NetworkModel}; is_per_unit::Bool=false)::Union{T,Tuple{T,MathBus,MathBranch}} where T <: MathSwitch
    math_obj = T(;
        index = id,
        f_bus = bus_lookup[eng_obj.f_bus],
        t_bus = bus_lookup[eng_obj.t_bus],
        f_connections = eng_obj.f_connections,
        t_connections = eng_obj.t_connections,
        status = Int(eng_obj.status),
        state = Int(eng_obj.state),
        dispatchable = Int(eng_obj.dispatchable),
        current_rating = eng_obj.cm_ub,
        thermal_rating = eng_obj.sm_ub,
    )

    if !(all(isapprox.(get(eng_obj, "rs", zeros(1, 1)), 0)) && all(isapprox.(get(eng_obj, "xs", zeros(1, 1)), 0)))
        bus_obj = MathBusObj(;
            index = length(eng.bus)+1,
            bus_i = length(eng.bus)+1,
            bus_type = math_obj.status == 0 ? 4 : 1,
            terminals = eng_obj.t_connections,
            grounded = fill(false, length(eng_obj.t_connections)),
            source_id = "$(typeof(T)).$(eng_obj.name)"
        )

        branch_obj = MathBranchObj(;
            index = length(eng.branch)+1,
            f_bus = math_obj.t_bus,
            t_bus = bus_obj.bus_i,
            f_connections = eng_obj.t_connections,
            t_connections = eng_obj.t_connections,
            br_r = eng_obj.rs .* eng_obj.length,
            br_x = eng_obj.xs .* eng_obj.length,
            angmin = fill(-10.0, length(eng_obj.t_connections)),
            angmax = fill( 10.0, length(eng_obj.t_connections)),
            c_rating_a = eng_obj.cm_ub,
            rate_a = eng_obj.sm_ub,
            br_status = Int(eng_obj.status)
        )

        math_obj.t_bus = bus_obj.bus_i

        return (math_obj, bus_obj, branch_obj)
    end

    return math_obj
end

function create_math_object(::Type{T}, eng_obj::EngTransformer, id::Int, bus_lookup::Dict{String,Int}, eng::EngineeringModel{NetworkModel}; is_per_unit::Bool=false, kron_reduced::Bool=false)::Union{T,Tuple{Vector{T},Vector{MathBus},Vector{MathBranch}}} where T <: MathTransformer
    transformers = MathTransformer[]
    buses = MathBus[]
    branches = MathBranch[]

    begin
        vnom = eng_obj.vm_nom * eng.settings.voltage_scale_factor
        snom = eng_obj.sm_nom * eng.settings.power_scale_factor

        nrw = length(eng_obj.bus)

        # calculate zbase in which the data is specified, and convert to SI
        zbase = (vnom.^2) ./ snom

        # x_sc is specified with respect to first winding
        x_sc = eng_obj.xsc .* zbase[1]

        # rs is specified with respect to each winding
        r_s = eng_obj.rw .* zbase

        g_sh =  (eng_obj.noloadloss*snom[1])/vnom[1]^2
        b_sh = -(eng_obj.cmag*snom[1])/vnom[1]^2

        # data is measured externally, but we now refer it to the internal side
        ratios = vnom/eng.settings.voltage_scale_factor
        x_sc = x_sc./ratios[1]^2
        r_s = r_s./ratios.^2
        g_sh = g_sh*ratios[1]^2
        b_sh = b_sh*ratios[1]^2

        # convert x_sc from list of upper triangle elements to an explicit dict
        y_sh = g_sh + im*b_sh
        z_sc = Dict([(key, im*x_sc[i]) for (i,key) in enumerate([(i,j) for i in 1:nrw for j in i+1:nrw])])

        dims = length(eng_obj.tm_set[1])

        (transformer_t_bus_w,transformer_t_branch_w) = build_loss_model(eng_obj.name, r_s, z_sc, y_sh,eng_obj.connections[1]; nphases=dims, status=Int(eng_obj.status), kron_reduced=kron_reduced)

        for w in 1:nrw
            # 2-WINDING TRANSFORMER
            # make virtual bus and mark it for reduction
            tm_nom = eng_obj.configurations[w]==DELTA ? eng_obj.vm_nom[w]*sqrt(3) : eng_obj.vm_nom[w]
            push!(transformers, T(;
                # name = "_virtual_transformer.$(eng_obj.name).$w",
                source_id = "_virtual_transformer.$(eng_obj.source_id).$w",
                f_bus = bus_lookup[eng_obj.bus[w]],
                t_bus = transformer_t_bus_w[w].bus_i,
                tm_nom = tm_nom,
                f_connections = eng_obj.connections[w],
                t_connections = kron_reduced ? eng_obj.connections[1] : collect(1:dims+1),
                configuration = eng_obj.configurations[w],
                polarity = eng_obj.polarity[w],
                tm_set = eng_obj.tm_set[w],
                tm_fix = eng_obj.tm_fix[w],
                sm_ub = eng_obj.sm_ub,
                cm_ub = eng_obj.cm_ub,
                tm_lb = eng_obj.tm_lb,
                tm_ub = eng_obj.tm_ub,
                tm_step = eng_obj.tm_step,
                status = Int(eng_obj.status),
                index = length(eng.transformer)+1
            ))

            # for prop in [["tm_lb", "tm_ub", "tm_step"]; pass_props]
            #     if haskey(eng_obj, prop)
            #         transformer_2wa_obj[prop] = eng_obj[prop][w]
            #     end
            # end

            # data_math["transformer"]["$(transformer_2wa_obj["index"])"] = transformer_2wa_obj

            # add regcontrol items to math model
            # if haskey(eng_obj,"controls") && !all(data_math["transformer"]["$(transformer_2wa_obj["index"])"]["tm_fix"])
            #     reg_obj = Dict{String,Any}(
            #         "vreg" => eng_obj.controls["vreg"][w],
            #         "band" => eng_obj.controls["band"][w],
            #         "ptratio" => eng_obj.controls["ptratio"][w],
            #         "ctprim" => eng_obj.controls["ctprim"][w],
            #         "r" => eng_obj.controls["r"][w],
            #         "x" => eng_obj.controls["x"][w],
            #     )
            #     data_math["transformer"]["$(transformer_2wa_obj["index"])"]["controls"] = reg_obj
            # end

            # if w==3 && eng_obj.polarity[w]==-1 # identify center-tapped transformer and mark all secondary-side nodes as triplex by adding va_start
            #     default_va = [0, -120, 120][eng_obj.connections[1][1]]
            #     data_math["bus"]["$(transformers[w].f_bus)"]["va_start"] = haskey(eng["bus"][eng_obj.bus[w]],"va_start") ? eng["bus"][eng_obj.bus[w]]["va_start"] : [default_va, (default_va+180)]
            #     idx = 0
            #     bus_ids = []
            #     t_bus = haskey(eng, "line") ? [data["t_bus"] for (_,data) in eng["line"] if data["f_bus"] == eng_obj.bus[w]] : []
            #     while length(t_bus)>0 || idx<length(bus_ids)
            #         for bus_idx in t_bus
            #             bus_id = data_math["bus_lookup"]["$bus_idx"]
            #             push!(bus_ids, bus_id)
            #             default_va = [0, -120, 120][eng_obj.connections[1][1]]
            #             data_math["bus"]["$bus_id"]["va_start"] = haskey(eng["bus"]["$bus_idx"],"va_start") ? eng["bus"]["$bus_idx"]["va_start"] : [default_va, (default_va+180)]
            #         end
            #         idx += 1
            #         t_bus = [data["t_bus"] for (_,data) in eng["line"] if data["f_bus"] == data_math["bus"]["$(bus_ids[idx])"]["name"]]
            #     end
            # end

            # push!(to_map, "transformer.$(transformer_2wa_obj["index"])")

        end
    end

    return (transformers, buses, branches)
end
