"""
    apply_kron_reduction!(data::Dict{String,<:Any}; kr_phases::Union{Vector{Int},Vector{String}}=[1,2,3], kr_neutral::Union{Int,String}=4)

Applies a Kron Reduction to the network, reducing out the `kr_neutral`, leaving only the `kr_phases`
"""
function apply_kron_reduction!(data::Dict{String,<:Any}; kr_phases::Union{Vector{Int},Vector{String}}=[1, 2, 3], kr_neutral::Union{Int,String}=4)
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_apply_kron_reduction!, data; apply_to_subnetworks=true, kr_phases=kr_phases, kr_neutral=kr_neutral)
end


"""
    _apply_kron_reduction!(data_eng::Dict{String,<:Any}; kr_phases::Union{Vector{Int},Vector{String}}=[1,2,3], kr_neutral::Union{Int,String}=4)

Applies a Kron Reduction to the network, reducing out the `kr_neutral`, leaving only the `kr_phases`
"""
function _apply_kron_reduction!(data_eng::Dict{String,<:Any}; kr_phases::Union{Vector{Int},Vector{String}}=[1, 2, 3], kr_neutral::Union{Int,String}=4)
    if !get(data_eng, "is_kron_reduced", false)
        if haskey(data_eng, "bus")
            for (id, eng_obj) in data_eng["bus"]
                filter = eng_obj["terminals"] .!= kr_neutral
                terminals_kr = eng_obj["terminals"][filter]

                @assert all(t in kr_phases for t in terminals_kr) "bus $id has terminals $(eng_obj["terminals"]), outside of $kr_phases, cannot be kron reduced"

                _apply_filter!(eng_obj, ["vm", "va", "vm_lb", "vm_ub", "vm_start", "va_start"], filter)
                eng_obj["terminals"] = terminals_kr

                gr_filter = eng_obj["grounded"] .!= kr_neutral
                _apply_filter!(eng_obj, ["grounded", "rg", "xg"], gr_filter)
            end
        end

        if haskey(data_eng, "line")
            for (_, eng_obj) in data_eng["line"]
                @assert all(eng_obj["f_connections"] .== eng_obj["t_connections"]) "Kron reduction is only supported if f_connections == t_connections"

                _apply_linecode!(eng_obj, data_eng)

                filter = _kron_reduce_branch!(eng_obj, ["rs", "xs"], ["g_fr", "b_fr", "g_to", "b_to"], eng_obj["f_connections"], kr_neutral)
                _apply_filter!(eng_obj, ["vad_lb", "vad_ub", "cm_ub", "sm_ub", "f_connections", "t_connections"], filter)
            end
        end

        if haskey(data_eng, "transformer")
            for (id, eng_obj) in data_eng["transformer"]
                _apply_xfmrcode!(eng_obj, data_eng)

                if haskey(eng_obj, "f_connections")
                    @assert all(eng_obj["f_connections"] .== eng_obj["t_connections"]) "Kron reduction is only supported if f_connections == t_connections"

                    filter = eng_obj["f_connections"] .!= kr_neutral
                    _apply_filter!(eng_obj, ["f_connections", "t_connections"], filter)
                else
                    for (w, connections) in enumerate(eng_obj["connections"])
                        filter = connections .!= kr_neutral
                        _apply_filter!(eng_obj, ["connections"], w, filter)
                    end
                end
            end
        end

        if haskey(data_eng, "switch")
            for (_, eng_obj) in data_eng["switch"]
                @assert all(eng_obj["f_connections"] .== eng_obj["t_connections"]) "Kron reduction is only supported if f_connections == t_connections"

                _apply_linecode!(eng_obj, data_eng)

                filter = _kron_reduce_branch!(eng_obj, Vector{String}([k for k in ["rs", "xs"] if haskey(eng_obj, k)]), String[], eng_obj["f_connections"], kr_neutral)
                _apply_filter!(eng_obj, ["vad_lb", "vad_ub", "cm_ub", "sm_ub", "f_connections", "t_connections"], filter)
            end
        end

        if haskey(data_eng, "shunt")
            for (_, eng_obj) in data_eng["shunt"]
                filter = _kron_reduce_branch!(eng_obj, String[], ["gs", "bs"], eng_obj["connections"], kr_neutral)
                _apply_filter!(eng_obj, ["connections"], filter)
            end
        end

        if haskey(data_eng, "load")
            for (_, eng_obj) in data_eng["load"]
                if eng_obj["configuration"] == WYE
                    @assert eng_obj["connections"][end] == kr_neutral "for wye-connected loads, to kron reduce the connections list should end with a neutral"

                    filter = eng_obj["connections"] .!= kr_neutral
                    _apply_filter!(eng_obj, ["connections"], filter)
                end
            end
        end

        if haskey(data_eng, "generator")
            for (_, eng_obj) in data_eng["generator"]
                if eng_obj["configuration"] == WYE
                    @assert eng_obj["connections"][end] == kr_neutral "for wye-connected generators, to kron reduce the connections list should end with a neutral"

                    filter = eng_obj["connections"] .!= kr_neutral
                    _apply_filter!(eng_obj, ["connections"], filter)
                end
            end
        end

        if haskey(data_eng, "solar")
            for (_, eng_obj) in data_eng["solar"]
                if eng_obj["configuration"] == WYE
                    @assert eng_obj["connections"][end] == kr_neutral "for wye-connected solar, to kron reduce the connections list should end with a neutral"

                    filter = eng_obj["connections"] .!= kr_neutral
                    _apply_filter!(eng_obj, ["connections"], filter)
                end
            end
        end

        if haskey(data_eng, "storage")
            for (_, eng_obj) in data_eng["storage"]
                filter = eng_obj["connections"] .!= kr_neutral
                _apply_filter!(eng_obj, ["connections"], filter)
            end
        end

        if haskey(data_eng, "voltage_source")
            for (_, eng_obj) in data_eng["voltage_source"]
                filter = eng_obj["connections"] .!= kr_neutral
                _apply_filter!(eng_obj, ["vm", "va", "vm_lb", "vm_ub", "rs", "xs", "connections"], filter)
            end
        end

        delete!(data_eng, "linecode")  # kron reduction moves linecode properties directly to lines
        delete!(data_eng, "xfmrcode")  # kron reduction moves xfmrcode properties directly to transformers

        find_conductor_ids!(data_eng)
        data_eng["is_kron_reduced"] = true
    end
end


"Return the Kron-reduction of the specified neutral conductors of a series impedance matrix."
function _kron_reduce_series_impedance(Z::Matrix, neutral_conductors::Vector{Int})
    # neutral conductor idxs
    N = neutral_conductors
    # phase conductor idxs (complement)
    P = setdiff(1:size(Z)[1], N)

    Zkr = Z[P, P] - Z[P, N] * inv(Z[N, N]) * Z[N, P]

    return Zkr
end


"Return the Kron-reduction of the specified neutral conductors of a shunt addmittance matrix."
function _kron_reduce_shunt_addmittance(Y::Matrix, neutral_conductors::Vector{Int})
    # phase conductor idxs (complement)
    P = setdiff(1:size(Y)[1], neutral_conductors)

    Ykr = Y[P, P]

    return Ykr
end


"Kron-reduce specified neutral conductors of a linecode."
function _kron_reduce_linecode!(l, neutral_conductors::Vector{Int})
    z_s = _kron_reduce_series_impedance(l["rs"] .+ im * l["xs"], neutral_conductors)
    y_fr = _kron_reduce_shunt_addmittance(l["g_fr"] .+ im * l["b_fr"], neutral_conductors)
    y_to = _kron_reduce_shunt_addmittance(l["g_to"] .+ im * l["b_to"], neutral_conductors)
    l["rs"] = real.(z_s)
    l["xs"] = imag.(z_s)
    l["g_fr"] = real.(y_fr)
    l["b_fr"] = imag.(y_fr)
    l["g_to"] = real.(y_to)
    l["b_to"] = imag.(y_to)
end


"""
    kron_reduce_implicit_neutrals!(data::Dict{String,Any})::Dict{String,Any}

Kron-reduce all (implied) neutral conductors of lines, switches and shunts, and remove any terminals which become unconnected.
A line or switch conductor is considered as a neutral conductor if it is connected between two neutral terminals.
A terminal is a neutral terminals if it is galvanically connected (i.e. through a line or switch)
to a grounded terminal, or the neutral conductor of a wye-connected component.
"""
function kron_reduce_implicit_neutrals!(data::Dict{String,Any})::Dict{String,Any}
    @assert iseng(data) "wrong data model type"

    apply_pmd!(_kron_reduce_implicit_neutrals!, data; apply_to_subnetworks=true)
end


"""
    _kron_reduce_implicit_neutrals!(data_eng::Dict{String,Any})::Dict{String,Any}

Kron-reduce all (implied) neutral conductors of lines, switches and shunts, and remove any terminals which become unconnected.
A line or switch conductor is considered as a neutral conductor if it is connected between two neutral terminals.
A terminal is a neutral terminals if it is galvanically connected (i.e. through a line or switch)
to a grounded terminal, or the neutral conductor of a wye-connected component.
"""
function _kron_reduce_implicit_neutrals!(data_eng::Dict{String,Any})::Dict{String,Any}
    # store the original linecode ids to detect naming clashes
    orig_lc_ids = [keys(data_eng["linecode"])...]

    # obtain all (implicit) neutral terminals
    nbts = _infer_neutral_terminals(data_eng)

    # Kron-reduce each line if eligible
    for (id, line) in get(data_eng, "line", Dict())
        doubly_grounded = [(line["f_bus"], t_fr) in nbts && (line["t_bus"], t_to) in nbts for (t_fr, t_to) in zip(line["f_connections"], line["t_connections"])]
        if any(doubly_grounded)
            keep = (!).(doubly_grounded)
            neutral_conductors = findall(doubly_grounded)
            _apply_filter!(line, ["f_connections", "t_connections", "cm_ub", "cm_ub_b", "cm_ub_c"], keep)
            if haskey(line, "linecode")
                suffix = "_kr_" * join(findall(doubly_grounded), ".")
                lc_orig_id = line["linecode"]
                lc_kr_id = lc_orig_id * suffix
                @assert !(lc_kr_id in orig_lc_ids) "Kron-reduced linecode naming clashes with original linecode names."
                line["linecode"] = lc_kr_id
                if !haskey(data_eng["linecode"], lc_kr_id)
                    data_eng["linecode"][lc_kr_id] = lc = deepcopy(data_eng["linecode"][lc_orig_id])
                    _kron_reduce_linecode!(lc, neutral_conductors)
                end
            end
            if haskey(line, "rs")
                _kron_reduce_linecode!(line, neutral_conductors)
            end
        end

        return data_eng
    end

    # Kron-reduce each shunt if eligible
    for (id, shunt) in get(data_eng, "shunt", Dict())
        cond_is_neutral = [(shunt["bus"], t) in nbts for t in shunt["connections"]]
        cond_keep = (!).(cond_is_neutral)
        _apply_filter!(shunt, ["connections"], cond_keep)
        Ys = _kron_reduce_shunt_addmittance(shunt["gs"] .+ im * shunt["bs"], findall(cond_is_neutral))
        shunt["gs"] = real.(Ys)
        shunt["bs"] = imag.(Ys)
    end

    # Kron-reduce each switch if eligible
    for (id, switch) in get(data_eng, "switch", Dict())
        doubly_grounded = [(switch["f_bus"], t_fr) in nbts && (switch["t_bus"], t_to) in nbts for (t_fr, t_to) in zip(switch["f_connections"], switch["t_connections"])]
        if any(doubly_grounded)
            keep = (!).(doubly_grounded)
            neutral_conductors = findall(doubly_grounded)
            _apply_filter!(switch, ["f_connections", "t_connections", "cm_ub", "cm_ub_b", "cm_ub_c"], keep)
            if haskey(switch, "linecode")
                suffix = "_kr_" * join(findall(doubly_grounded), ".")
                lc_orig_id = switch["linecode"]
                lc_kr_id = lc_orig_id * suffix
                @assert !(lc_kr_id in orig_lc_ids) "Kron-reduced linecode naming clashes with original linecode names."
                line["linecode"] = lc_kr_id
                if !haskey(data_eng["linecode"], lc_kr_id)
                    data_eng["linecode"][lc_kr_id] = _kron_reduce_linecode(data_eng["linecode"][lc_orig_id], neutral_conductors)
                end
            end
            if haskey(switch, "rs")
                Zs = _kron_reduce_series_impedance(switch["rs"] .+ im * switch["xs"], findall(doubly_grounded))
                switch["rs"] = real.(Zs)
                switch["xs"] = imag.(Zs)
            end
        end
    end

    # remove unconnected terminals (likely to be caused by earlier Kron-reductions)
    remove_unconnected_terminals!(data_eng)

    # ground remaining neutral terminals
    remaining_bts = [(b, t) for (b, bus) in data_eng["bus"] for t in bus["terminals"]]
    for (b, t) in intersect(nbts, remaining_bts)
        bus = data_eng["bus"][b]
        perfectly_grounded = bus["grounded"][iszero.(bus["rg"] .+ im * bus["xg"])]
        if !(t in perfectly_grounded)
            push!(bus["grounded"], t)
            push!(bus["rg"], 0.0)
            push!(bus["xg"], 0.0)
        end
    end

    # update 'conductor_ids' property
    find_conductor_ids!(data_eng)

    return data_eng
end
