"items that are mapped one-to-one from engineering to math models"
const _1to1_maps = Dict{String,Vector{String}}(
    "bus" => ["vm", "va", "vm_start", "va_start", "terminals", "phases", "neutral", "vm_pn_lb", "vm_pn_ub", "vm_pp_lb", "vm_pp_ub", "vm_ng_ub", "dss", "vuf_ub", "vm_pair_lb", "vm_pair_ub"],
    "line" => ["f_connections", "t_connections", "dss"],
    "transformer" => ["f_connections", "t_connections", "dss"],
    "switch" => ["status", "f_connections", "t_connections", "dss"],
    "shunt" => ["status", "dispatchable", "gs", "bs", "connections", "dss"],
    "load" => ["model", "configuration", "connections", "dispatchable", "status", "dss"],
    "generator" => ["configuration", "connections", "dss"],
    "solar" => ["configuration", "connections", "dss"],
    "storage" => ["status", "energy", "configuration", "connections", "dss"],
    "voltage_source" => ["configuration", "connections", "dss"],
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
        kron_reduce::Bool=true,
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

If `kron_reduce==true`, [`apply_kron_reduction!`](@ref apply_kron_reduction!) and
[`apply_phase_projection_delta!`](@ref apply_phase_projection_delta!) will be applied to
the network data.

## Phase projection

If `phase_project==true`,  [`apply_phase_projection!`](@ref apply_phase_projection!) will be
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
ENGINEERING model, a property called `z` was added to `switch` objects, and a property at
the root level of the dictionary was added called `max_switch_actions`, the user could
pass the following dictionary to `eng2math_passthrough`:

    Dict{String,Vector{String}}(
        "switch" => String["z"],
        "root" => String["max_switch_actions"],
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
    kron_reduce::Bool=true,
    phase_project::Bool=false,
    multinetwork::Bool=false,
    global_keys::Set{String}=Set{String}(),
    eng2math_passthrough::Dict{String,<:Vector{<:String}}=Dict{String,Vector{String}}(),
    eng2math_extensions::Vector{<:Function}=Function[],
    make_pu::Bool=true,
    make_pu_extensions::Vector{<:Function}=Function[],
    correct_network_data::Bool=true,
    )::Dict{String,Any}

    current_data_model = get(data, "data_model", MATHEMATICAL)

    if iseng(data)
        if multinetwork && !ismultinetwork(data)
            data = make_multinetwork(data; global_keys=global_keys)
        end

        data_math = _map_eng2math(
            data;
            kron_reduce=kron_reduce,
            phase_project=phase_project,
            eng2math_extensions=eng2math_extensions,
            eng2math_passthrough=eng2math_passthrough,
            global_keys=global_keys,
        )
        correct_network_data && correct_network_data!(data_math; make_pu=make_pu, make_pu_extensions=make_pu_extensions)

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
    kron_reduce::Bool=true,
    phase_project::Bool=false,
    eng2math_extensions::Vector{<:Function}=Function[],
    eng2math_passthrough::Dict{String,Vector{String}}=Dict{String,Vector{String}}(),
    global_keys::Set{String}=Set{String}(),
    )::Dict{String,Any}

    @assert iseng(data_eng)

    _data_eng = deepcopy(data_eng)

    if kron_reduce
        apply_kron_reduction!(_data_eng)
        apply_phase_projection_delta!(_data_eng)
    end

    if phase_project
        apply_phase_projection!(_data_eng)
    end

    if ismultinetwork(data_eng)
        data_math = Dict{String,Any}(
            "name" => get(_data_eng, "name", ""),
            "data_model" => MATHEMATICAL,
            "nw" => Dict{String,Any}(
                n => Dict{String,Any}(
                    "per_unit" => get(_data_eng, "per_unit", false),
                    "is_projected" => get(nw, "is_projected", false),
                    "is_kron_reduced" => get(nw, "is_kron_reduced", false),
                    "settings" => deepcopy(nw["settings"]),
                    "time_elapsed" => get(nw, "time_elapsed", 1.0),
                ) for (n,nw) in _data_eng["nw"]
            ),
            "multinetwork" => ismultinetwork(data_eng),
            [k => data_eng[k] for k in global_keys if haskey(data_eng, k)]...
        )
    else
        data_math = Dict{String,Any}(
            "name" => get(_data_eng, "name", ""),
            "per_unit" => get(_data_eng, "per_unit", false),
            "data_model" => MATHEMATICAL,
            "is_projected" => get(_data_eng, "is_projected", false),
            "is_kron_reduced" => get(_data_eng, "is_kron_reduced", false),
            "settings" => deepcopy(_data_eng["settings"]),
            "time_elapsed" => get(_data_eng, "time_elapsed", 1.0),
        )
    end

    apply_pmd!(_map_eng2math_nw!, data_math, _data_eng; eng2math_passthrough=eng2math_passthrough, eng2math_extensions=eng2math_extensions)

    if ismultinetwork(data_eng)
        _collect_nw_maps!(data_math)
        _collect_nw_bus_lookups!(data_math)
    end

    return data_math
end


"""
"""
function _collect_nw_maps!(data_math::Dict{String,<:Any})
    @assert ismultinetwork(data_math)
    @assert ismath(data_math)

    data_math["map"] = Dict{String,Vector{Dict{String,Any}}}()
    for (n,nw) in data_math["nw"]
        data_math["map"][n] = pop!(nw, "map")
    end
end


"""
"""
function _collect_nw_bus_lookups!(data_math::Dict{String,<:Any})
    @assert ismultinetwork(data_math)
    @assert ismath(data_math)

    data_math["bus_lookup"] = Dict{String,Dict{String,Any}}()
    for (n,nw) in data_math["nw"]
        data_math["bus_lookup"][n] = pop!(nw, "bus_lookup")
    end
end


"""
"""
function _map_eng2math_nw!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; eng2math_passthrough::Dict{String,Vector{String}}=Dict{String,Vector{String}}(), eng2math_extensions::Vector{<:Function}=Function[])
    data_math["map"] = Vector{Dict{String,Any}}([
        Dict{String,Any}("unmap_function" => "_map_math2eng_root!")
    ])

    _init_base_components!(data_math)

    for property in get(eng2math_passthrough, "root", String[])
        if haskey(data_eng, property)
            data_math[property] = deepcopy(data_eng[property])
        end
    end

    for type in pmd_eng_asset_types
        getfield(PowerModelsDistribution, Symbol("_map_eng2math_$(type)!"))(data_math, data_eng; pass_props=get(eng2math_passthrough, type, String[]))
    end

    # Custom eng2math transformation functions
    for eng2math_func! in eng2math_extensions
        eng2math_func!(data_math, data_eng)
    end

    # post fix
    if !get(data_math, "is_kron_reduced", false)
        #TODO fix this in place / throw error instead? IEEE8500 leads to switches
        # with 3x3 R matrices but only 1 phase
        #NOTE: Don't do this when kron-reducing, it will undo the padding
        _slice_branches!(data_math)
    end

    find_conductor_ids!(data_math)
    _map_conductor_ids!(data_math)

    _map_settings_vbases_default!(data_math)
end


function _map_settings_vbases_default!(data_math::Dict{String,<:Any})
    vbases_default = Dict{String,Real}()
    for (bus,vbase) in get(data_math["settings"], "vbases_default", Dict())
        vbases_default["$(data_math["bus_lookup"][bus])"] = vbase
    end

    data_math["settings"]["vbases_default"] = vbases_default
end


"converts engineering bus components into mathematical bus components"
function _map_eng2math_bus!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    for (name, eng_obj) in get(data_eng, "bus", Dict{String,Any}())
        terminals = eng_obj["terminals"]

        math_obj = _init_math_obj("bus", name, eng_obj, length(data_math["bus"])+1; pass_props=pass_props)

        math_obj["bus_i"] = math_obj["index"]
        math_obj["bus_type"] = eng_obj["status"] == DISABLED ? 4 : 1
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
                "bs" => imag.(sh_y),
            )
            push!(to_sh, "shunt.$sh_index")
        end

        if haskey(eng_obj, "vm")
            math_obj["vm"] = eng_obj["vm"]
        end
        if haskey(eng_obj, "va")
            math_obj["va"] = eng_obj["va"]
        end

        math_obj["vmin"], math_obj["vmax"] = _get_tight_absolute_voltage_magnitude_bounds(eng_obj)
        math_obj["vm_pair_lb"], math_obj["vm_pair_ub"] = _get_tight_pairwise_voltage_magnitude_bounds(eng_obj)
        _add_implicit_absolute_bounds!(math_obj, terminals)

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
                math_obj[t_key] = haskey(eng_obj, f_key) ? eng_obj[f_key] : fill(Inf, nphases)
        end

        math_obj["br_status"] = eng_obj["status"] == DISABLED ? 0 : 1

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
                "sm_ub" => get(eng_obj, "sm_ub", Inf),
                "cm_ub" => get(eng_obj, "cm_ub", Inf),
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

            @info "---- NAME --- : $(name)"
            @info "RS: $(r_s)"
            @info "XSC: $(x_sc)"
            @info "GSH: $(g_sh)"
            @info "BSH: $(b_sh)"

            # convert x_sc from list of upper triangle elements to an explicit dict
            y_sh = g_sh + im*b_sh
            z_sc = Dict([(key, im*x_sc[i]) for (i,key) in enumerate([(i,j) for i in 1:nrw for j in i+1:nrw])])

            @info "$(y_sh)"
            @info "$(z_sc)"

            dims = length(eng_obj["tm_set"][1])
            transformer_t_bus_w = _build_loss_model!(data_math, name, to_map, r_s, z_sc, y_sh,eng_obj["connections"][1]; nphases=dims, status=Int(eng_obj["status"] == ENABLED))

            for w in 1:nrw
                # 2-WINDING TRANSFORMER
                # make virtual bus and mark it for reduction
                tm_nom = eng_obj["configuration"][w]==DELTA ? eng_obj["vm_nom"][w]*sqrt(3) : eng_obj["vm_nom"][w]

                @info "TM_NOM: $(tm_nom)"

                transformer_2wa_obj = Dict{String,Any}(
                    "name"          => "_virtual_transformer.$name.$w",
                    "source_id"     => "_virtual_transformer.$(eng_obj["source_id"]).$w",
                    "f_bus"         => data_math["bus_lookup"][eng_obj["bus"][w]],
                    "t_bus"         => transformer_t_bus_w[w],
                    "tm_nom"        => tm_nom,
                    "f_connections" => eng_obj["connections"][w],
                    "t_connections" => get(data_math, "is_kron_reduced", false) ? eng_obj["connections"][1] : collect(1:dims+1),
                    "configuration" => eng_obj["configuration"][w],
                    "polarity"      => eng_obj["polarity"][w],
                    "tm_set"        => eng_obj["tm_set"][w],
                    "tm_fix"        => eng_obj["tm_fix"][w],
                    "sm_ub"         => get(eng_obj, "sm_ub", Inf),
                    "cm_ub"         => get(eng_obj, "cm_ub", Inf),
                    "status"        => eng_obj["status"] == DISABLED ? 0 : 1,
                    "index"         => length(data_math["transformer"])+1
                )

                for prop in [["tm_lb", "tm_ub", "tm_step"]; pass_props]
                    if haskey(eng_obj, prop)
                        transformer_2wa_obj[prop] = eng_obj[prop][w]
                    end
                end

                data_math["transformer"]["$(transformer_2wa_obj["index"])"] = transformer_2wa_obj

                # add regcontrol items to math model
                if haskey(eng_obj,"controls") && !all(data_math["transformer"]["$(transformer_2wa_obj["index"])"]["tm_fix"])
                    reg_obj = Dict{String,Any}(
                        "vreg" => eng_obj["controls"]["vreg"][w],
                        "band" => eng_obj["controls"]["band"][w],
                        "ptratio" => eng_obj["controls"]["ptratio"][w],
                        "ctprim" => eng_obj["controls"]["ctprim"][w],
                        "r" => eng_obj["controls"]["r"][w],
                        "x" => eng_obj["controls"]["x"][w],
                    )
                    data_math["transformer"]["$(transformer_2wa_obj["index"])"]["controls"] = reg_obj
                end
                if w==3 && eng_obj["polarity"][w]==-1 # identify center-tapped transformer and mark all secondary-side nodes as triplex by adding va_start
                    default_va = [0, -120, 120][eng_obj["connections"][1][1]]
                    data_math["bus"]["$(transformer_2wa_obj["f_bus"])"]["va_start"] = haskey(data_eng["bus"][eng_obj["bus"][w]],"va_start") ? data_eng["bus"][eng_obj["bus"][w]]["va_start"] : [default_va, (default_va+180)]
                    idx = 0
                    bus_ids = []
                    t_bus = haskey(data_eng, "line") ? [data["t_bus"] for (_,data) in data_eng["line"] if data["f_bus"] == eng_obj["bus"][w]] : []
                    while length(t_bus)>0 || idx<length(bus_ids)
                        for bus_idx in t_bus
                            bus_id = data_math["bus_lookup"]["$bus_idx"]
                            push!(bus_ids, bus_id)
                            default_va = [0, -120, 120][eng_obj["connections"][1][1]]
                            data_math["bus"]["$bus_id"]["va_start"] = haskey(data_eng["bus"]["$bus_idx"],"va_start") ? data_eng["bus"]["$bus_idx"]["va_start"] : [default_va, (default_va+180)]
                        end
                        idx += 1
                        t_bus = [data["t_bus"] for (_,data) in data_eng["line"] if data["f_bus"] == data_math["bus"]["$(bus_ids[idx])"]["name"]]
                    end
                end

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
        math_obj["status"] = eng_obj["status"] == DISABLED ? 0 : 1

        math_obj["state"] = Int(get(eng_obj, "state", CLOSED))
        math_obj["dispatchable"] = Int(get(eng_obj, "dispatchable", YES))

        # OPF bounds
        for (f_key, t_key) in [("cm_ub", "current_rating"), ("cm_ub_b", "c_rating_b"), ("cm_ub_c", "c_rating_c"),
            ("sm_ub", "thermal_rating"), ("sm_ub_b", "rate_b"), ("sm_ub_c", "rate_c")]
            math_obj[t_key] = haskey(eng_obj, f_key) ? eng_obj[f_key] : fill(Inf, nphases)
        end

        map_to = "switch.$(math_obj["index"])"

        if haskey(eng_obj, "linecode")
            _apply_linecode!(eng_obj, data_eng)
        end

        if !(all(isapprox.(get(eng_obj, "rs", zeros(1, 1)), 0)) && all(isapprox.(get(eng_obj, "xs", zeros(1, 1)), 0)))
            # build virtual bus

            f_bus = deepcopy(data_math["bus"]["$(math_obj["f_bus"])"])
            t_bus = deepcopy(data_math["bus"]["$(math_obj["t_bus"])"])

            N = length(eng_obj["t_connections"])
            bus_obj = Dict{String,Any}(
                "name" => "_virtual_bus.switch.$name",
                "bus_i" => length(data_math["bus"])+1,
                "bus_type" => eng_obj["status"] == DISABLED ? 4 : 1,
                "terminals" => eng_obj["t_connections"],  # connected to the switch on the to-side
                "grounded" => fill(false, N),  # connected to the switch on the to-side
                "vmin" => fill(0.0, N),
                "vmax" => fill(Inf, N),
                "vm_pair_lb" => Tuple{Any,Any,Real}[],
                "vm_pair_ub" => Tuple{Any,Any,Real}[],
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
                "g_fr" => _admittance_conversion(data_eng, eng_obj, "g_fr"),
                "g_to" => _admittance_conversion(data_eng, eng_obj, "g_to"),
                "b_fr" => _admittance_conversion(data_eng, eng_obj, "b_fr"),
                "b_to" => _admittance_conversion(data_eng, eng_obj, "b_to"),
                "angmin" => fill(-10.0, nphases),
                "angmax" => fill( 10.0, nphases),
                "c_rating_a" => fill(Inf, nphases),
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

        # add capcontrol items to math model
        if haskey(eng_obj,"controls")
            math_obj["controls"] = deepcopy(eng_obj["controls"])
            dss_obj_type = split(math_obj["controls"]["element"], "."; limit=2)[1]
            if dss_obj_type == "line"
                elem_id = filter(x->x.second["source_id"] == math_obj["controls"]["element"], data_math["branch"])
                if !isempty(elem_id)
                    elem_id = first(elem_id).first

                    math_obj["controls"]["element"] = Dict{String,Any}(
                        "type" => "branch",
                        "index" => data_math["branch"][elem_id]["index"],
                        "f_bus" => data_math["branch"][elem_id]["f_bus"],
                        "t_bus" => data_math["branch"][elem_id]["t_bus"]
                    )
                else
                    elem_id = first(filter(x->x.second["source_id"] == replace(math_obj["controls"]["element"], "line."=>"switch."), data_math["switch"])).first
                    math_obj["controls"]["element"] = Dict{String,Any}(
                        "type" => "switch",
                        "index" => data_math["switch"][elem_id]["index"],
                        "f_bus" => data_math["switch"][elem_id]["f_bus"],
                        "t_bus" => data_math["switch"][elem_id]["t_bus"]
                    )
                end
            elseif dss_obj_type == "capacitor"
                elem_id = first(filter(x->x.second["source_id"] == replace(math_obj["controls"]["element"], "capacitor"=>"shunt"), data_math["shunt"])).first
                math_obj["controls"]["element"] = Dict{String,Any}(
                    "type" => "shunt",
                    "index" => data_math["shunt"][elem_id]["index"],
                    "f_bus" => data_math["shunt"][elem_id]["shunt_bus"],
                    "t_bus" => data_math["shunt"][elem_id]["shunt_bus"]
                )
            else
                elem_id = first(filter(x->x.second["source_id"] == math_obj["controls"]["element"], data_math["transformer"])).first
                math_obj["controls"]["element"] = Dict{String,Any}(
                    "type" => "transformer",
                    "index" => data_math["transformer"][elem_id]["index"],
                    "f_bus" => data_math["transformer"][elem_id]["f_bus"],
                    "t_bus" => data_math["transformer"][elem_id]["t_bus"]
                )
            end
        end

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
        if eng_obj["model"]==ZIP
            to_map = String[]
            for (idx,l) in enumerate([IMPEDANCE, CURRENT, POWER])
                math_obj = Dict{String,Any}(
                    "model" => l,
                    "connections" => eng_obj["connections"],
                    "configuration" => eng_obj["configuration"],
                    "name" => "$(name)_$l",
                    "status" => eng_obj["status"] == ENABLED ? 1 : 0,
                    "qd" => eng_obj["qd_nom"]*eng_obj["zipv"][3+idx],
                    "vnom_kv" => eng_obj["vm_nom"],
                    "source_id" => eng_obj["source_id"],
                    "load_bus" => data_math["bus_lookup"][eng_obj["bus"]],
                    "dispatchable" => eng_obj["dispatchable"] == NO ? 0 : 1,
                    "index" => length(data_math["load"])+1,
                    "pd" => eng_obj["pd_nom"]*eng_obj["zipv"][idx],
                )

                data_math["load"]["$(math_obj["index"])"] = math_obj

                push!(to_map, "load.$(math_obj["index"])")
            end
            push!(data_math["map"], Dict{String,Any}(
                "from" => name,
                "to" => to_map,
                "unmap_function" => "_map_math2eng_load!",
            ))
        else
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
end


"converts engineering generators into mathematical generators"
function _map_eng2math_generator!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any}; pass_props::Vector{String}=String[])
    for (name, eng_obj) in get(data_eng, "generator", Dict{String,Any}())
        math_obj = _init_math_obj("generator", name, eng_obj, length(data_math["gen"])+1; pass_props=pass_props)

        connections = eng_obj["connections"]

        math_obj["gen_bus"] = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["gen_status"] = status = Int(eng_obj["status"])
        math_obj["control_mode"] = control_mode = Int(get(eng_obj, "control_mode", FREQUENCYDROOP))
        math_obj["pmax"] = get(eng_obj, "pg_ub", fill(Inf, length(connections)))

        bus_type = data_math["bus"]["$(math_obj["gen_bus"])"]["bus_type"]
        data_math["bus"]["$(math_obj["gen_bus"])"]["bus_type"] = _compute_bus_type(bus_type, status, control_mode)
        if control_mode == Int(ISOCHRONOUS) && status == 1
            data_math["bus"]["$(math_obj["gen_bus"])"]["vm"] = eng_obj["vg"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["vmax"] = eng_obj["vg"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["vmin"] = eng_obj["vg"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["va"] = [0.0, -120, 120, zeros(length(data_math["bus"]["$(math_obj["gen_bus"])"]) - 3)...][data_math["bus"]["$(math_obj["gen_bus"])"]["terminals"]]
        end

        for (f_key, t_key) in [("qg_lb", "qmin"), ("qg_ub", "qmax"), ("pg_lb", "pmin")]
            if haskey(eng_obj, f_key)
                math_obj[t_key] = eng_obj[f_key]
            elseif f_key[end-1:end]=="ub"
                math_obj[t_key] = fill(Inf, length(math_obj["pmax"]))
            else
                math_obj[t_key] = fill(-Inf, length(math_obj["pmax"]))
            end
        end

        math_obj["pg"] = get(eng_obj, "pg", fill(0.0, length(connections)))
        math_obj["qg"] = get(eng_obj, "qg", fill(0.0, length(connections)))

        _add_gen_cost_model!(math_obj, eng_obj)

        math_obj["configuration"] = get(eng_obj, "configuration", WYE)

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

        connections = eng_obj["connections"]

        math_obj["gen_bus"] = data_math["bus_lookup"][eng_obj["bus"]]
        math_obj["gen_status"] = status = Int(eng_obj["status"])

        math_obj["control_mode"] = control_mode = Int(get(eng_obj, "control_mode", FREQUENCYDROOP))
        bus_type = data_math["bus"]["$(math_obj["gen_bus"])"]["bus_type"]
        data_math["bus"]["$(math_obj["gen_bus"])"]["bus_type"] = _compute_bus_type(bus_type, status, control_mode)
        if control_mode == Int(ISOCHRONOUS) && status == 1
            data_math["bus"]["$(math_obj["gen_bus"])"]["vm"] = eng_obj["vg"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["vmax"] = eng_obj["vg"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["vmin"] = eng_obj["vg"]
            data_math["bus"]["$(math_obj["gen_bus"])"]["va"] = [0.0, -120, 120, zeros(length(data_math["bus"]["$(math_obj["gen_bus"])"]) - 3)...][data_math["bus"]["$(math_obj["gen_bus"])"]["terminals"]]
            data_math["bus"]["$(math_obj["gen_bus"])"]["bus_type"] = 3
        end

        for (fr_k, to_k) in [("vg", "vg")]
            if haskey(eng_obj, fr_k)
                math_obj[to_k] = eng_obj[fr_k]
            end
        end

        N = eng_obj["configuration"]==DELTA && length(eng_obj["connections"])==1 ? 1 : _infer_int_dim_unit(eng_obj, false) # if solar is delta-connected to triplex node, N can be equal to 1
        for (fr_k, to_k, def) in [("pg_lb", "pmin", -Inf), ("pg_ub", "pmax", Inf), ("qg_lb", "qmin", -Inf), ("qg_ub", "qmax", Inf)]
            math_obj[to_k] = haskey(eng_obj, fr_k) ? eng_obj[fr_k] : fill(def, N)
        end

        math_obj["pg"] = get(eng_obj, "pg", fill(0.0, N))
        math_obj["qg"] = get(eng_obj, "qg", fill(0.0, N))

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
        math_obj["configuration"] = get(eng_obj, "configuration", WYE)

        math_obj["energy"] = eng_obj["energy"]
        math_obj["energy_rating"] = eng_obj["energy_ub"]
        math_obj["charge_rating"] = eng_obj["charge_ub"]
        math_obj["discharge_rating"] = eng_obj["discharge_ub"]
        math_obj["charge_efficiency"] = eng_obj["charge_efficiency"] / 100.0
        math_obj["discharge_efficiency"] = eng_obj["discharge_efficiency"] / 100.0
        math_obj["thermal_rating"] = get(eng_obj, "sm_ub", Inf)
        math_obj["qmin"] = eng_obj["qs_lb"]
        math_obj["qmax"] = eng_obj["qs_ub"]
        math_obj["r"] = eng_obj["rs"]
        math_obj["x"] = eng_obj["xs"]
        math_obj["p_loss"] = eng_obj["pex"]
        math_obj["q_loss"] = eng_obj["qex"]

        math_obj["ps"] = get(eng_obj, "ps", 0.0)
        math_obj["qs"] = get(eng_obj, "qs", 0.0)

        math_obj["control_mode"] = control_mode = Int(get(eng_obj, "control_mode", FREQUENCYDROOP))
        bus_type = data_math["bus"]["$(math_obj["storage_bus"])"]["bus_type"]
        data_math["bus"]["$(math_obj["storage_bus"])"]["bus_type"] = _compute_bus_type(bus_type, math_obj["status"], control_mode)
        if control_mode == Int(ISOCHRONOUS) && math_obj["status"] == 1
            data_math["bus"]["$(math_obj["storage_bus"])"]["va"] = [0.0, -120, 120, zeros(length(data_math["bus"]["$(math_obj["storage_bus"])"]) - 3)...][data_math["bus"]["$(math_obj["storage_bus"])"]["terminals"]]
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
        math_obj["gen_status"] = status = Int(eng_obj["status"])
        math_obj["pg"] = get(eng_obj, "pg", fill(0.0, nphases))
        math_obj["qg"] = get(eng_obj, "qg", fill(0.0, nphases))
        math_obj["vg"] = eng_obj["vm"]
        math_obj["pmin"] = get(eng_obj, "pg_lb", fill(-Inf, nphases))
        math_obj["pmax"] = get(eng_obj, "pg_ub", fill( Inf, nphases))
        math_obj["qmin"] = get(eng_obj, "qg_lb", fill(-Inf, nphases))
        math_obj["qmax"] = get(eng_obj, "qg_ub", fill( Inf, nphases))
        math_obj["connections"] = eng_obj["connections"]
        math_obj["configuration"] = get(eng_obj, "configuration", WYE)
        math_obj["control_mode"] = control_mode = Int(get(eng_obj, "control_mode", ISOCHRONOUS))
        math_obj["source_id"] = "voltage_source.$name"

        _add_gen_cost_model!(math_obj, eng_obj)

        map_to = "gen.$(math_obj["index"])"

        if !all(isapprox.(get(eng_obj, "rs", zeros(1, 1)), 0)) && !all(isapprox.(get(eng_obj, "xs", zeros(1, 1)), 0))
            f_bus = deepcopy(data_math["bus"]["$(math_obj["gen_bus"])"])

            bus_obj = Dict{String,Any}(
                "bus_i" => length(data_math["bus"])+1,
                "index" => length(data_math["bus"])+1,
                "terminals" => eng_obj["connections"],
                "grounded" => [f_bus["grounded"][findfirst(f_bus["terminals"].==t)] for t in eng_obj["connections"]],
                "name" => "_virtual_bus.voltage_source.$name",
                "bus_type" => status == 0 ? 4 : control_mode == Int(ISOCHRONOUS) ? 3 : 2,
                "vm" => deepcopy(eng_obj["vm"]),
                "va" => deepcopy(eng_obj["va"]),
                "vmin" => deepcopy(get(eng_obj, "vm_lb", control_mode == Int(ISOCHRONOUS) ? eng_obj["vm"] : fill(0.0, nphases))),
                "vmax" => deepcopy(get(eng_obj, "vm_ub", control_mode == Int(ISOCHRONOUS) ? eng_obj["vm"] : fill(Inf, nphases))),
                "vm_pair_lb" => deepcopy(get(eng_obj, "vm_pair_lb", Tuple{Any,Any,Real}[])),
                "vm_pair_ub" => deepcopy(get(eng_obj, "vm_pair_ub", Tuple{Any,Any,Real}[])),
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
                "c_rating_a" => fill(Inf, nconductors),
                "br_status" => status,
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
            vm_lb = control_mode == Int(ISOCHRONOUS) ? eng_obj["vm"] : get(eng_obj, "vm_lb", fill(0.0, nphases))
            vm_ub = control_mode == Int(ISOCHRONOUS) ? eng_obj["vm"] : get(eng_obj, "vm_ub", fill(Inf, nphases))

            data_math["bus"]["$gen_bus"]["vmin"] = [vm_lb..., [0.0 for n in 1:(nconductors-nphases)]...]
            data_math["bus"]["$gen_bus"]["vmax"] = [vm_ub..., [Inf for n in 1:(nconductors-nphases)]...]
            data_math["bus"]["$gen_bus"]["vm"] = [eng_obj["vm"]..., [0.0 for n in 1:(nconductors-nphases)]...]
            data_math["bus"]["$gen_bus"]["va"] = [eng_obj["va"]..., [0.0 for n in 1:(nconductors-nphases)]...]

            bus_type = data_math["bus"]["$gen_bus"]["bus_type"]
            data_math["bus"]["$gen_bus"]["bus_type"] = _compute_bus_type(bus_type, status, control_mode)
        end

        data_math["gen"]["$(math_obj["index"])"] = math_obj

        push!(data_math["map"], Dict{String,Any}(
            "from" => name,
            "to" => map_to,
            "unmap_function" => "_map_math2eng_voltage_source!",
        ))
    end
end
