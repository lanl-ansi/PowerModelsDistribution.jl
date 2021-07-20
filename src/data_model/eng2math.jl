import LinearAlgebra: diagm

"items that are mapped one-to-one from engineering to math models"
const _1to1_maps = Dict{String,Vector{String}}(
    "bus" => ["vm", "va", "vm_start", "va_start", "terminals", "phases", "neutral", "vm_pn_lb", "vm_pn_ub", "vm_pp_lb", "vm_pp_ub", "vm_ng_ub", "dss", "vuf_ub", "vm_pair_lb", "vm_pair_ub"],
    "line" => ["f_connections", "t_connections", "dss"],
    "transformer" => ["f_connections", "t_connections", "dss"],
    "switch" => ["status", "f_connections", "t_connections", "dss"],
    "shunt" => ["status", "dispatchable", "gs", "bs", "connections", "dss"],
    "load" => ["model", "configuration", "connections", "dispatchable", "status", "dss"],
    "generator" => ["pg", "qg", "vg", "configuration", "connections", "dss"],
    "solar" => ["pg", "qg", "configuration", "connections", "dss"],
    "storage" => ["status", "energy", "ps", "qs", "connections", "dss"],
    "voltage_source" => ["dss"],
)

"list of nodal type elements in the engineering model"
const _eng_node_elements = String[
    "load", "shunt", "generator", "solar", "storage", "voltage_source"
]

"list of edge type elements in the engineering model"
const _eng_edge_elements = String[
    "line", "switch", "transformer"
]

"list of all eng asset types"
const pmd_eng_asset_types = String[
    "bus", _eng_edge_elements..., _eng_node_elements...
]

"list of nodal type elements in the engineering model"
const _math_node_elements = String[
    "load", "shunt", "gen", "storage"
]

"list of edge type elements in the engineering model"
const _math_edge_elements = String[
    "branch", "switch", "transformer"
]

"list of math asset types that are dispatchable"
const _math_dispatchable_elements = String[
    "load", "shunt", "switch"
]

"list of all math asset types"
const pmd_math_asset_types = String[
    "bus", _math_node_elements..., _math_edge_elements...
]


"""
    transform_data_model(
        data::Dict{String,<:Any};
        kron_reduced::Bool=true,
        multinetwork::Bool=false,
        global_keys::Set{String}=Set{String}(),
        eng2math_passthrough::Dict{String,<:Vector{<:String}}=Dict{String,Vector{String}}(),
        eng2math_extensions::Vector{<:Function}=Function[],
        make_pu::Bool=true,
        make_pu_extensions::Vector{<:Function}=Function[],
    )::Dict{String,Any}

Transforms a data model model between ENGINEERING (high-level) and MATHEMATICAL (low-level)
[`DataModel`](@ref DataModel).

# Notes

## Kron reduction

If `kron_reduced==true`, [`apply_kron_reduction!`](@ref apply_kron_reduction!) will be
applied to the network data.

## Multinetwork transformations

If `multinetwork==true`, the data model will be transformed into a multinetwork (e.g.,
time series) data structure using [`make_multinetwork`](@ref make_multinetwork) before
being transformed into a MATHEMATICAL [`DataModel`](@ref DataModel).

`global_keys::Set{String}` can be used to add custom top-level items to the multinetwork
data structure, and will only be used in the context where `multinetwork==true`, and
ignored otherwise.

## Custom eng2math transformations

To add custom transformations between ENGINEERING and MATHEMATICAL data models,
`eng2math_extensions::Vector{<:Function}` can be utilized to pass user-created functions,
which are expected to have the signature

    eng2math_func!(data_math::Dict{String,Any}, data_eng::Dict{String,Any})

where data_math and data_eng equivalent to single subnetworks in a multinetwork data structure,
or a non-multinetwork data structure.

These functions are run after all built-in eng2math transformations have been performed.

### Mapping back to ENGINEERING

See [`transform_solution`](@ref transform_solution)

### Passthrough properties

To more simply pass through some properties in the built-in eng2math transformations,
`eng2math_passthrough::Dict{String,Vector{String}}` can be used. For example, if in the
ENGINEERING model, a property called `z` was added to `switch` objects, the user could
pass the following dictionary to `eng2math_passthrough`:

    Dict{String,Vector{String}}(
        "switch" => String["z"],
    )

This will result in `z` showing up on the `switch` object in the MATHEMATICAL model.
Passthrough properties will always land on the __primary__ conversion object in the
MATHEMATICAL model if that object gets converted to multiple object types, e.g.,
`voltage_source` with internal impedance will result in `gen`, `bus`, and `branch`
objects in the MATHEMATICAL model, but passthrough properties will only land on `gen`.

## Custom per-unit transformations

To add additional per-unit transformations, a user can supply custom functions to
`make_pu_extensions::Vector{<:Function}`, which will only be used if `make_pu==true`.

See [`make_per_unit!`](@ref make_per_unit!) for further explanation.
"""
function transform_data_model(
    data::Dict{String,<:Any};
    kron_reduced::Bool=true,
    multinetwork::Bool=false,
    global_keys::Set{String}=Set{String}(),
    eng2math_passthrough::Dict{String,<:Vector{<:String}}=Dict{String,Vector{String}}(),
    eng2math_extensions::Vector{<:Function}=Function[],
    make_pu::Bool=true,
    make_pu_extensions::Vector{<:Function}=Function[],
    )::Dict{String,Any}

    current_data_model = get(data, "data_model", MATHEMATICAL)

    if iseng(data)
        if multinetwork && !ismultinetwork(data)
            data = make_multinetwork(data; global_keys=global_keys)
        end

        data_math = _map_eng2math(
            data;
            kron_reduced=kron_reduced,
            eng2math_extensions=eng2math_extensions,
            eng2math_passthrough=eng2math_passthrough
        )
        correct_network_data!(data_math; make_pu=make_pu, make_pu_extensions=make_pu_extensions)

        return data_math
    elseif ismath(data)
        @info "A MATHEMATICAL data model cannot be converted back to an ENGINEERING data model, irreversible transformations have already been made"
        return data
    else
        @info "Data model '$current_data_model' is not recognized, no model type transformation performed"
        return data
    end
end


"base function for converting engineering model to mathematical model"
function _map_eng2math(
    data_eng::Dict{String,<:Any};
    kron_reduced::Bool=true,
    eng2math_extensions::Vector{<:Function}=Function[],
    eng2math_passthrough::Dict{String,Vector{String}}=Dict{String,Vector{String}}(),
    global_keys::Set{String}=Set{String}(),
    )::Dict{String,Any}

    @assert iseng(data_eng)

    _data_eng = deepcopy(data_eng)

    data_math = Dict{String,Any}(
        "name" => get(_data_eng, "name", ""),
        "per_unit" => get(_data_eng, "per_unit", false),
        "data_model" => MATHEMATICAL,
        "nw" => Dict{String,Any}(),
        "multinetwork" => ismultinetwork(data_eng)
    )

    if ismultinetwork(data_eng)
        nw_data_eng = _data_eng["nw"]
    else
        nw_data_eng = Dict("0" => _data_eng)
    end

    for (n, nw_eng) in nw_data_eng
        # TODO remove kron reduction from eng2math (breaking)
        if kron_reduced && !get(nw_eng, "is_kron_reduced", false)
            apply_kron_reduction!(nw_eng)
        end

        if !get(data_eng, "is_projected", false)
            apply_phase_projection_delta!(nw_eng)
        end

        nw_math = Dict{String,Any}(
            "is_projected" => get(nw_eng, "is_projected", false),
            "is_kron_reduced" => get(nw_eng, "is_kron_reduced", false),
            "settings" => deepcopy(nw_eng["settings"]),
        )

        if haskey(nw_eng, "time_elapsed")
            nw_math["time_elapsed"] = nw_eng["time_elapsed"]
        end

        #TODO the PM tests break for branches which are not of the size indicated by conductors;
        # for now, set to 1 to prevent this from breaking when not kron-reduced

        nw_math["map"] = Vector{Dict{String,Any}}([
            Dict{String,Any}("unmap_function" => "_map_math2eng_root!")
        ])

        _init_base_components!(nw_math)

        for type in pmd_eng_asset_types
            getfield(PowerModelsDistribution, Symbol("_map_eng2math_$(type)!"))(nw_math, nw_eng; pass_props=get(eng2math_passthrough, type, String[]))
        end

        # Custom eng2math transformation functions
        for eng2math_func! in eng2math_extensions
            eng2math_func!(nw_math, nw_eng)
        end

        # post fix
        if !get(nw_math, "is_kron_reduced", false)
            #TODO fix this in place / throw error instead? IEEE8500 leads to switches
            # with 3x3 R matrices but only 1 phase
            #NOTE: Don't do this when kron-reducing, it will undo the padding
            _slice_branches!(nw_math)
        end

        find_conductor_ids!(nw_math)
        _map_conductor_ids!(nw_math)

        data_math["nw"][n] = nw_math
    end

    for k in union(_pmd_math_global_keys, global_keys)
        for (n,nw) in data_math["nw"]
            if haskey(nw, k)
                data_math[k] = pop!(nw, k)
            end
        end
    end

    if !ismultinetwork(data_eng)
        merge!(data_math, pop!(data_math["nw"], "0"))
        delete!(data_math, "nw")
        delete!(data_math, "multinetwork")
    end

    return data_math
end


"converts engineering bus components into mathematical bus components"
function _map_eng2math_bus!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    for (name, eng_obj) in get(data_eng, "bus", Dict{String,Any}())
        terminals = eng_obj["terminals"]

        math_obj = _init_math_obj("bus", name, eng_obj, length(data_math["bus"])+1; pass_props=pass_props)

        math_obj["bus_i"] = math_obj["index"]
        math_obj["bus_type"] = _bus_type_conversion(data_eng, eng_obj, "status")
        math_obj["source_id"] = "bus.$name"

        # take care of grounding; convert to shunt if lossy
        grounded_perfect, shunts = _convert_grounding(eng_obj["terminals"], eng_obj["grounded"], eng_obj["rg"], eng_obj["xg"])

        math_obj["grounded"] = grounded_perfect
        to_sh = []
        for (sh_connections, sh_y) in shunts
            sh_index = length(data_math["shunt"]) + 1
            data_math["shunt"]["$sh_index"] = Dict(
                "index" => sh_index,
                "shunt_bus" => math_obj["bus_i"],
                "connections" => sh_connections,
                "gs" => real.(sh_y),
                "bs" => real.(sh_y),
            )
            push!(to_sh, "shunt.$sh_index")
        end

        if haskey(eng_obj, "vm")
            math_obj["vm"] = eng_obj["vm"]
        end
        if haskey(eng_obj, "va")
            math_obj["va"] = eng_obj["va"]
        end

        math_obj["vmin"] = get(eng_obj, "vm_lb", fill(0.0, length(terminals)))
        math_obj["vmax"] = get(eng_obj, "vm_ub", fill(Inf, length(terminals)))

        data_math["bus"]["$(math_obj["index"])"] = math_obj

        if !haskey(data_math, "bus_lookup")
            data_math["bus_lookup"] = Dict{Any,Int}()
        end

        data_math["bus_lookup"][name] = math_obj["index"]

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "bus.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_bus!",
        ))
    end
end


"converts engineering lines into mathematical branches"
function _map_eng2math_line!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    for (name, eng_obj) in get(data_eng, "line", Dict{Any,Dict{String,Any}}())
        _apply_linecode!(eng_obj, data_eng)

        math_obj = _init_math_obj("line", name, eng_obj, length(data_math["branch"])+1; pass_props=pass_props)

        nphases = size(eng_obj["rs"])[1]

        math_obj["f_bus"] = data_math["bus_lookup"][eng_obj["f_bus"]]
        math_obj["t_bus"] = data_math["bus_lookup"][eng_obj["t_bus"]]

        math_obj["br_r"] = _impedance_conversion(data_eng, eng_obj, "rs")
        math_obj["br_x"] = _impedance_conversion(data_eng, eng_obj, "xs")

        math_obj["g_fr"] = _admittance_conversion(data_eng, eng_obj, "g_fr")
        math_obj["g_to"] = _admittance_conversion(data_eng, eng_obj, "g_to")

        math_obj["b_fr"] = _admittance_conversion(data_eng, eng_obj, "b_fr")
        math_obj["b_to"] = _admittance_conversion(data_eng, eng_obj, "b_to")

        math_obj["angmin"] = get(eng_obj, "vad_lb", fill(-60.0, nphases))
        math_obj["angmax"] = get(eng_obj, "vad_ub", fill( 60.0, nphases))

        for (f_key, t_key) in [("cm_ub", "c_rating_a"), ("cm_ub_b", "c_rating_b"), ("cm_ub_c", "c_rating_c"),
            ("sm_ub", "rate_a"), ("sm_ub_b", "rate_b"), ("sm_ub_c", "rate_c")]
            if haskey(eng_obj, f_key)
                math_obj[t_key] = eng_obj[f_key]
            end
        end

        math_obj["br_status"] = Int(eng_obj["status"])

        data_math["branch"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "branch.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_line!",
        ))
    end
end


"converts engineering n-winding transformers into mathematical ideal 2-winding lossless transformer branches and impedance branches to represent the loss model"
function _map_eng2math_transformer!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    for (name, eng_obj) in get(data_eng, "transformer", Dict{Any,Dict{String,Any}}())
        # Build map first, so we can update it as we decompose the transformer
        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => String[],
            "unmap_function" => "_map_math2eng_transformer!",
        ))

        to_map = data_math["map"][end]["to"]

        _apply_xfmrcode!(eng_obj, data_eng)

        if haskey(eng_obj, "f_bus") && haskey(eng_obj, "t_bus")
            @assert all(haskey(eng_obj, k) for k in ["f_bus", "t_bus", "f_connections", "t_connections"]) "Incomplete definition of AL2W tranformer $name, aborting eng2math conversion"

            nphases = length(eng_obj["f_connections"])

            math_obj = Dict{String,Any}(
                "name" => name,
                "source_id" => eng_obj["source_id"],
                "f_bus" => data_math["bus_lookup"][eng_obj["f_bus"]],
                "t_bus" => data_math["bus_lookup"][eng_obj["t_bus"]],
                "f_connections" => eng_obj["f_connections"],
                "t_connections" => eng_obj["t_connections"],
                "configuration" => get(eng_obj, "configuration", WYE),
                "tm_nom" => get(eng_obj, "tm_nom", 1.0),
                "tm_set" => get(eng_obj, "tm_set", fill(1.0, nphases)),
                "tm_fix" => get(eng_obj, "tm_fix", fill(true, nphases)),
                "polarity" => get(eng_obj, "polarity", -1),
                "status" => Int(get(eng_obj, "status", ENABLED)),
                "index" => length(data_math["transformer"])+1
            )

            for k in [["tm_lb", "tm_ub"]; pass_props]
                if haskey(eng_obj, k)
                    math_obj[k] = eng_obj[k]
                end
            end

            data_math["transformer"]["$(math_obj["index"])"] = math_obj

            push!(to_map, "transformer.$(math_obj["index"])")
        else
            vnom = eng_obj["vm_nom"] * data_eng["settings"]["voltage_scale_factor"]
            snom = eng_obj["sm_nom"] * data_eng["settings"]["power_scale_factor"]

            nrw = length(eng_obj["bus"])

            # calculate zbase in which the data is specified, and convert to SI
            zbase = (vnom.^2) ./ snom

            # x_sc is specified with respect to first winding
            x_sc = eng_obj["xsc"] .* zbase[1]

            # rs is specified with respect to each winding
            r_s = eng_obj["rw"] .* zbase

            g_sh =  (eng_obj["noloadloss"]*snom[1])/vnom[1]^2
            b_sh = -(eng_obj["cmag"]*snom[1])/vnom[1]^2

            # data is measured externally, but we now refer it to the internal side
            ratios = vnom/data_eng["settings"]["voltage_scale_factor"]
            x_sc = x_sc./ratios[1]^2
            r_s = r_s./ratios.^2
            g_sh = g_sh*ratios[1]^2
            b_sh = b_sh*ratios[1]^2

            # convert x_sc from list of upper triangle elements to an explicit dict
            y_sh = g_sh + im*b_sh
            z_sc = Dict([(key, im*x_sc[i]) for (i,key) in enumerate([(i,j) for i in 1:nrw for j in i+1:nrw])])

            dims = length(eng_obj["tm_set"][1])
            transformer_t_bus_w = _build_loss_model!(data_math, name, to_map, r_s, z_sc, y_sh; nphases=dims)

            for w in 1:nrw
                # 2-WINDING TRANSFORMER
                # make virtual bus and mark it for reduction
                tm_nom = eng_obj["configuration"][w]==DELTA ? eng_obj["vm_nom"][w]*sqrt(3) : eng_obj["vm_nom"][w]
                transformer_2wa_obj = Dict{String,Any}(
                    "name"          => "_virtual_transformer.$name.$w",
                    "source_id"     => "_virtual_transformer.$(eng_obj["source_id"]).$w",
                    "f_bus"         => data_math["bus_lookup"][eng_obj["bus"][w]],
                    "t_bus"         => transformer_t_bus_w[w],
                    "tm_nom"        => tm_nom,
                    "f_connections" => eng_obj["connections"][w],
                    "t_connections" => get(data_math, "is_kron_reduced", false) ? collect(1:dims) : collect(1:dims+1),
                    "configuration" => eng_obj["configuration"][w],
                    "polarity"      => eng_obj["polarity"][w],
                    "tm_set"        => eng_obj["tm_set"][w],
                    "tm_fix"        => eng_obj["tm_fix"][w],
                    "status"        => Int(get(eng_obj, "status", ENABLED)),
                    "index"         => length(data_math["transformer"])+1
                )

                for prop in [["tm_lb", "tm_ub", "tm_step"]; pass_props]
                    if haskey(eng_obj, prop)
                        transformer_2wa_obj[prop] = eng_obj[prop][w]
                    end
                end

                data_math["transformer"]["$(transformer_2wa_obj["index"])"] = transformer_2wa_obj

                push!(to_map, "transformer.$(transformer_2wa_obj["index"])")
            end
        end
    end
end


"converts engineering switches into mathematical switches and (if neeed) impedance branches to represent loss model"
function _map_eng2math_switch!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    # TODO enable real switches (right now only using vitual lines)
    for (name, eng_obj) in get(data_eng, "switch", Dict{Any,Dict{String,Any}}())
        nphases = length(eng_obj["f_connections"])

        math_obj = _init_math_obj("switch", name, eng_obj, length(data_math["switch"])+1; pass_props=pass_props)

        math_obj["f_bus"] = data_math["bus_lookup"][eng_obj["f_bus"]]
        math_obj["t_bus"] = data_math["bus_lookup"][eng_obj["t_bus"]]

        math_obj["state"] = Int(get(eng_obj, "state", CLOSED))
        math_obj["dispatchable"] = Int(get(eng_obj, "dispatchable", YES))

        # OPF bounds
        for (f_key, t_key) in [("cm_ub", "c_rating_a"), ("cm_ub_b", "c_rating_b"), ("cm_ub_c", "c_rating_c"),
            ("sm_ub", "rate_a"), ("sm_ub_b", "rate_b"), ("sm_ub_c", "rate_c")]
            if haskey(eng_obj, f_key)
                math_obj[t_key] = eng_obj[f_key]
            end
        end

        map_to = "switch.$(math_obj["index"])"

        if haskey(eng_obj, "linecode")
            _apply_linecode!(eng_obj, data_eng)
        end

        if !(all(isapprox.(get(eng_obj, "rs", zeros(1, 1)), 0)) && all(isapprox.(get(eng_obj, "xs", zeros(1, 1)), 0)))
            # build virtual bus

            f_bus = deepcopy(data_math["bus"]["$(math_obj["f_bus"])"])
            t_bus = deepcopy(data_math["bus"]["$(math_obj["t_bus"])"])

            bus_obj = Dict{String,Any}(
                "name" => "_virtual_bus.switch.$name",
                "bus_i" => length(data_math["bus"])+1,
                "bus_type" => get(eng_obj, "status", ENABLED) == DISABLED ? 4 : 1,
                "terminals" => t_bus["terminals"],  # connected to the switch on the to-side
                "grounded" => t_bus["grounded"],  # connected to the switch on the to-side
                "vmin" => t_bus["vmin"],
                "vmax" => t_bus["vmax"],
                "source_id" => "switch.$name",
                "index" => length(data_math["bus"])+1,
            )

            math_obj["t_bus"] = bus_obj["bus_i"]
            data_math["bus"]["$(bus_obj["index"])"] = bus_obj

            branch_obj = _init_math_obj("line", name, eng_obj, length(data_math["branch"])+1)

            _branch_obj = Dict{String,Any}(
                "name" => "_virtual_branch.switch.$name",
                "source_id" => "switch.$name",
                "f_bus" => bus_obj["bus_i"],
                "t_bus" => data_math["bus_lookup"][eng_obj["t_bus"]],
                "f_connections" => eng_obj["t_connections"],  # the virtual branch connects to the switch on the to-side
                "t_connections" => eng_obj["t_connections"],  # should be identical to the switch's to-side connections
                "br_r" => _impedance_conversion(data_eng, eng_obj, "rs"),
                "br_x" => _impedance_conversion(data_eng, eng_obj, "xs"),
                "g_fr" => zeros(nphases, nphases),
                "g_to" => zeros(nphases, nphases),
                "b_fr" => zeros(nphases, nphases),
                "b_to" => zeros(nphases, nphases),
                "angmin" => fill(-10.0, nphases),
                "angmax" => fill( 10.0, nphases),
                "br_status" => eng_obj["status"] == DISABLED ? 0 : 1,
            )

            merge!(branch_obj, _branch_obj)

            data_math["branch"]["$(branch_obj["index"])"] = branch_obj

            map_to = [map_to, "bus.$(bus_obj["index"])", "branch.$(branch_obj["index"])"]
        end

        data_math["switch"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => map_to,
            "unmap_function" => "_map_math2eng_switch!",
        ))
    end
end


"converts engineering generic shunt components into mathematical shunt components"
function _map_eng2math_shunt!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    for (name, eng_obj) in get(data_eng, "shunt", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("shunt", name, eng_obj, length(data_math["shunt"])+1; pass_props=pass_props)

        # TODO change to new capacitor shunt calc logic
        math_obj["shunt_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["gs"] = get(eng_obj, "gs", zeros(size(eng_obj["bs"])))

        data_math["shunt"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "shunt.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_shunt!",
        ))
    end
end


"converts engineering load components into mathematical load components"
function _map_eng2math_load!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    for (name, eng_obj) in get(data_eng, "load", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("load", name, eng_obj, length(data_math["load"])+1; pass_props=pass_props)

        connections = eng_obj["connections"]

        math_obj["load_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["pd"] = eng_obj["pd_nom"]
        math_obj["qd"] = eng_obj["qd_nom"]

        math_obj["vnom_kv"] = eng_obj["vm_nom"]

        data_math["load"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "load.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_load!",
        ))
   end
end


"converts engineering generators into mathematical generators"
function _map_eng2math_generator!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    for (name, eng_obj) in get(data_eng, "generator", Dict{String,Any}())
        math_obj = _init_math_obj("generator", name, eng_obj, length(data_math["gen"])+1; pass_props=pass_props)

        connections = eng_obj["connections"]

        math_obj["gen_bus"] = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["gen_status"] = Int(eng_obj["status"])
        math_obj["control_mode"] = get(eng_obj, "control_mode", FREQUENCYDROOP)
        math_obj["pmax"] = get(eng_obj, "pg_ub", fill(Inf, length(connections)))

        math_obj["control_mode"] = control_mode = get(eng_obj, "control_mode", FREQUENCYDROOP)
        data_math["bus"]["$(math_obj["gen_bus"])"]["bus_type"] = control_mode == ISOCHRONOUS ? 3 : 2
        if control_mode == ISOCHRONOUS
            data_math["bus"]["$(math_obj["gen_bus"])"]["vm"] = eng_obj["vg"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["vmax"] = eng_obj["vg"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["vmin"] = eng_obj["vg"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["va"] = [0.0, -120, 120, zeros(length(data_math["bus"]["$(math_obj["gen_bus"])"]) - 3)...]
        end

        for (f_key, t_key) in [("qg_lb", "qmin"), ("qg_ub", "qmax"), ("pg_lb", "pmin")]
            if haskey(eng_obj, f_key)
                math_obj[t_key] = eng_obj[f_key]
            end
        end

        _add_gen_cost_model!(math_obj, eng_obj)

        math_obj["configuration"] = get(eng_obj, "configuration", WYE)

        # if PV generator mode convert attached bus to PV bus
        if math_obj["control_mode"] == ISOCHRONOUS
            data_math["bus"]["$(data_math["bus_lookup"][eng_obj["bus"]])"]["bus_type"] = 2
        end

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "gen.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_generator!",
        ))
    end
end


"converts engineering solar components into mathematical generators"
function _map_eng2math_solar!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    for (name, eng_obj) in get(data_eng, "solar", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("solar", name, eng_obj, length(data_math["gen"])+1; pass_props=pass_props)

        math_obj["gen_bus"] = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["gen_status"] = Int(eng_obj["status"])

        math_obj["control_mode"] = control_mode = get(eng_obj, "control_mode", FREQUENCYDROOP)
        data_math["bus"]["$(math_obj["gen_bus"])"]["bus_type"] = control_mode == ISOCHRONOUS ? 3 : 2
        if control_mode == ISOCHRONOUS
            data_math["bus"]["$(math_obj["gen_bus"])"]["vm"] = eng_obj["vg"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["vmax"] = eng_obj["vg"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["vmin"] = eng_obj["vg"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["va"] = [0.0, -120, 120, zeros(length(data_math["bus"]["$(math_obj["gen_bus"])"]) - 3)...]
        end

        for (fr_k, to_k) in [("vg", "vg"), ("pg_lb", "pmin"), ("pg_ub", "pmax"), ("qg_lb", "qmin"), ("qg_ub", "qmax")]
            if haskey(eng_obj, fr_k)
                math_obj[to_k] = eng_obj[fr_k]
            end
        end

        _add_gen_cost_model!(math_obj, eng_obj)

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "gen.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_solar!",
        ))
    end
end


"converts engineering storage into mathematical storage"
function _map_eng2math_storage!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    for (name, eng_obj) in get(data_eng, "storage", Dict{Any,Dict{String,Any}}())
        math_obj = _init_math_obj("storage", name, eng_obj, length(data_math["storage"])+1; pass_props=pass_props)

        math_obj["storage_bus"] = data_math["bus_lookup"][eng_obj["bus"]]

        math_obj["energy"] = eng_obj["energy"]
        math_obj["energy_rating"] = eng_obj["energy_ub"]
        math_obj["charge_rating"] = eng_obj["charge_ub"]
        math_obj["discharge_rating"] = eng_obj["discharge_ub"]
        math_obj["charge_efficiency"] = eng_obj["charge_efficiency"] / 100.0
        math_obj["discharge_efficiency"] = eng_obj["discharge_efficiency"] / 100.0
        math_obj["thermal_rating"] = eng_obj["cm_ub"]
        math_obj["qmin"] = eng_obj["qs_lb"]
        math_obj["qmax"] =  eng_obj["qs_ub"]
        math_obj["r"] = eng_obj["rs"]
        math_obj["x"] = eng_obj["xs"]
        math_obj["p_loss"] = eng_obj["pex"]
        math_obj["q_loss"] = eng_obj["qex"]

        math_obj["ps"] = get(eng_obj, "ps", zeros(size(eng_obj["cm_ub"])))
        math_obj["qs"] = get(eng_obj, "qs", zeros(size(eng_obj["cm_ub"])))

        math_obj["control_mode"] = control_mode = get(eng_obj, "control_mode", FREQUENCYDROOP)
        data_math["bus"]["$(math_obj["storage_bus"])"]["bus_type"] = control_mode == ISOCHRONOUS ? 3 : 2
        if control_mode == ISOCHRONOUS
            data_math["bus"]["$(math_obj["storage_bus"])"]["va"] = [0.0, -120, 120, zeros(length(data_math["bus"]["$(math_obj["gen_bus"])"]) - 3)...]
        end

        data_math["storage"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => "storage.$(math_obj["index"])",
            "unmap_function" => "_map_math2eng_storage!",
        ))
    end
end


"converts engineering voltage sources into mathematical generators and (if needed) impedance branches to represent the loss model"
function _map_eng2math_voltage_source!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    for (name, eng_obj) in get(data_eng, "voltage_source", Dict{String,Any}())
        nconductors = length(eng_obj["connections"])
        nphases = get(eng_obj, "configuration", WYE) == WYE && !get(data_eng, "is_kron_reduced", false) ? nconductors - 1 : nconductors

        math_obj = _init_math_obj("voltage_source", name, eng_obj, length(data_math["gen"])+1; pass_props=pass_props)

        math_obj["name"] = "_virtual_gen.voltage_source.$name"
        math_obj["gen_bus"] = gen_bus = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["connections"] = eng_obj["connections"]
        math_obj["gen_status"] = Int(eng_obj["status"])
        math_obj["pg"] = fill(0.0, nphases)
        math_obj["qg"] = fill(0.0, nphases)
        math_obj["vg"] = eng_obj["vm"]
        math_obj["pmin"] = get(eng_obj, "pg_lb", fill(-Inf, nphases))
        math_obj["pmax"] = get(eng_obj, "pg_ub", fill( Inf, nphases))
        math_obj["qmin"] = get(eng_obj, "qg_lb", fill(-Inf, nphases))
        math_obj["qmax"] = get(eng_obj, "qg_ub", fill( Inf, nphases))
        math_obj["connections"] = eng_obj["connections"]
        math_obj["configuration"] = get(eng_obj, "configuration", WYE)
        math_obj["control_mode"] = get(eng_obj, "control_mode", ISOCHRONOUS)
        math_obj["source_id"] = "voltage_source.$name"

        _add_gen_cost_model!(math_obj, eng_obj)

        map_to = "gen.$(math_obj["index"])"

        if !all(isapprox.(get(eng_obj, "rs", zeros(1, 1)), 0)) && !all(isapprox.(get(eng_obj, "xs", zeros(1, 1)), 0))
            f_bus = deepcopy(data_math["bus"]["$(math_obj["gen_bus"])"])

            bus_obj = Dict{String,Any}(
                "bus_i" => length(data_math["bus"])+1,
                "index" => length(data_math["bus"])+1,
                "terminals" => f_bus["terminals"],
                "grounded" => f_bus["grounded"],
                "name" => "_virtual_bus.voltage_source.$name",
                "bus_type" => math_obj["control_mode"] == ISOCHRONOUS ? 3 : 2,
                "vm" => deepcopy(eng_obj["vm"]),
                "va" => deepcopy(eng_obj["va"]),
                "vmin" => deepcopy(get(eng_obj, "vm_lb", math_obj["control_mode"] == ISOCHRONOUS ? eng_obj["vm"] : fill(0.0, nphases))),
                "vmax" => deepcopy(get(eng_obj, "vm_ub", math_obj["control_mode"] == ISOCHRONOUS ? eng_obj["vm"] : fill(Inf, nphases))),
                "source_id" => "voltage_source.$name",
            )
            for (i,t) in enumerate(eng_obj["connections"])
                if data_math["bus"]["$(data_math["bus_lookup"][eng_obj["bus"]])"]["grounded"][i]
                    bus_obj["vm"][i] = 0
                    bus_obj["vmin"][i] = 0
                    bus_obj["vmax"][i] = Inf
                end
            end

            math_obj["gen_bus"] = gen_bus = bus_obj["bus_i"]

            data_math["bus"]["$(bus_obj["index"])"] = bus_obj

            branch_obj = Dict{String,Any}(
                "name" => "_virtual_branch.voltage_source.$name",
                "source_id" => "voltage_source.$name",
                "f_bus" => bus_obj["bus_i"],
                "t_bus" => data_math["bus_lookup"][eng_obj["bus"]],
                "f_connections" => eng_obj["connections"],
                "t_connections" => eng_obj["connections"],
                "angmin" => fill(-10.0, nconductors),
                "angmax" => fill( 10.0, nconductors),
                "br_status" => 1,
                "br_r" => _impedance_conversion(data_eng, eng_obj, "rs"),
                "br_x" => _impedance_conversion(data_eng, eng_obj, "xs"),
                "g_fr" => zeros(nconductors, nconductors),
                "g_to" => zeros(nconductors, nconductors),
                "b_fr" => zeros(nconductors, nconductors),
                "b_to" => zeros(nconductors, nconductors),
                "index" => length(data_math["branch"])+1
            )


            data_math["branch"]["$(branch_obj["index"])"] = branch_obj

            map_to = [map_to, "bus.$(bus_obj["index"])", "branch.$(branch_obj["index"])"]
        else
            vm_lb = math_obj["control_mode"] == ISOCHRONOUS ? eng_obj["vm"] : get(eng_obj, "vm_lb", fill(0.0, nphases))
            vm_ub = math_obj["control_mode"] == ISOCHRONOUS ? eng_obj["vm"] : get(eng_obj, "vm_ub", fill(Inf, nphases))

            data_math["bus"]["$gen_bus"]["vmin"] = [vm_lb..., [0.0 for n in 1:(nconductors-nphases)]...]
            data_math["bus"]["$gen_bus"]["vmax"] = [vm_ub..., [Inf for n in 1:(nconductors-nphases)]...]
            data_math["bus"]["$gen_bus"]["vm"] = [eng_obj["vm"]..., [0.0 for n in 1:(nconductors-nphases)]...]
            data_math["bus"]["$gen_bus"]["va"] = [eng_obj["va"]..., [0.0 for n in 1:(nconductors-nphases)]...]
            data_math["bus"]["$gen_bus"]["bus_type"] = math_obj["control_mode"] == ISOCHRONOUS ? 3 : 2
        end

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => map_to,
            "unmap_function" => "_map_math2eng_voltage_source!",
        ))
    end
end
