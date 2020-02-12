#TODO
# Can buses in a voltage zone have different terminals?
# Add current/power bounds to data model


function copy_kwargs_to_dict_if_present!(dict, kwargs, args)
    for arg in args
        if haskey(kwargs, arg)
            dict[string(arg)] = kwargs[arg]
        end
    end
end

function add_kwarg!(dict, kwargs, name, default)
    if haskey(kwargs, name)
        dict[string(name)] = kwargs[name]
    else
        dict[string(name)] = default
    end
end

function component_dict_from_list!(list)
    dict = Dict{String, Any}()
    for comp_dict in list
        dict[comp_dict["id"]] = comp_dict
    end
    return dict
end

DTYPES = Dict{Symbol, Any}()
CHECKS = Dict{Symbol, Any}()

function check_dtypes(dict, dtypes, comp_type, id)
    for key in keys(dict)
        symb = Symbol(key)
        if haskey(dtypes, symb)
            @assert(isa(dict[key], dtypes[symb]), "$comp_type $id: the property $key should be a $(dtypes[symb]), not a $(typeof(dict[key])).")
        else
            @assert(false, "$comp_type $id: the property $key is unknown.")
        end
    end
end

function check_data_model(data)
    for component in [:bus, :linecode, :line, :load, :generator, :transformer_nph3w_lossy, :transformer_2w_ideal, :capacitor, :shunt]
        if haskey(data, string(component))
            for (id, comp_dict) in data[string(component)]
                if haskey(DTYPES, component)
                    check_dtypes(comp_dict, DTYPES[component], component, id)
                end
                if haskey(CHECKS, component)
                    CHECKS[component](data, comp_dict)
                end
            end
        end
    end
end


function create_data_model(;kwargs...)
    data_model = Dict{String, Any}()
    data_model["quantity_scalars"] = Dict{String, Any}()
    add_kwarg!(data_model, kwargs, :v_var_scalar, 1E3)
    return data_model
end

# COMPONENTS
#*#################

DTYPES[:ev] = Dict()
DTYPES[:storage] = Dict()
DTYPES[:pv] = Dict()
DTYPES[:wind] = Dict()
DTYPES[:switch] = Dict()
DTYPES[:shunt] = Dict()
DTYPES[:autotransformer] = Dict()
DTYPES[:synchronous_generator] = Dict()
DTYPES[:zip_load] = Dict()
DTYPES[:grounding] = Dict(
    :bus => String,
    :rg => Real,
    :xg => Real,
)
DTYPES[:synchronous_generator] = Dict()
DTYPES[:boundary] = Dict()
DTYPES[:meter] = Dict()



# Linecode

DTYPES[:linecode] = Dict(
    :id => AbstractString,
    :n_conductors => Int,
    :rs => Array{<:Real, 2},
    :xs => Array{<:Real, 2},
    :g_fr => Array{<:Real, 2},
    :g_to => Array{<:Real, 2},
    :b_fr => Array{<:Real, 2},
    :b_to => Array{<:Real, 2},
)


CHECKS[:linecode] = function check_line(data, linecode)
end


function create_linecode(id, n_conductors; kwargs...)
    linecode = Dict{String,Any}()
    linecode["id"] = id
    linecode["n_conductors"] = n_conductors
    add_kwarg!(linecode, kwargs, :rs, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(linecode, kwargs, :xs, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(linecode, kwargs, :g_fr, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(linecode, kwargs, :b_fr, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(linecode, kwargs, :g_to, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(linecode, kwargs, :b_to, fill(0.0, n_conductors, n_conductors))
    return linecode
end


# line

DTYPES[:line] = Dict(
    :id => AbstractString,
    :status => Int,
    :f_bus => String,
    :t_bus => String,
    :n_conductors => Int,
    :f_terminals => Vector{<:Int},
    :t_terminals => Vector{<:Int},
    :linecode => Int,
    :length => Real,
    :rs =>Matrix{<:Real},
    :xs =>Matrix{<:Real},
    :g_fr => Matrix{<:Real},
    :b_fr =>Matrix{<:Real},
    :g_to =>Matrix{<:Real},
    :b_to =>Matrix{<:Real},
    :b_to =>Matrix{<:Real},
    :c_rating =>Vector{<:Real},
)


CHECKS[:line] = function check_line(data, line)
    i = line["id"]

    if haskey(line, "linecode")
        # line is defined with a linecode
        @assert(haskey(line, "length"), "line $i: a line defined through a linecode, should have a length property.")

        linecode_id = line["linecode"]
        @assert(haskey(data, "linecode") && haskey(data["linecode"], "$linecode_id"), "line $i: the linecode $linecode_id is not defined.")
        linecode = data["linecode"]["$linecode_id"]

        for key in ["n_conductors", "rs", "xs", "g_fr", "g_to", "b_fr", "b_to"]
            @assert(!haskey(line, key), "line $i: a line with a linecode, should not specify $key; this is already done by the linecode.")
        end

        N = linecode["n_conductors"]
        @assert(length(line["f_terminals"])==N, "line $i: the number of terminals should match the number of conductors in the linecode.")
        @assert(length(line["t_terminals"])==N, "line $i: the number of terminals should match the number of conductors in the linecode.")
    else
        # normal line
        @assert(!haskey(line, "length"), "line $i: length only makes sense for linees defined through linecodes.")
        for key in ["n_conductors", "rs", "xs", "g_fr", "g_to", "b_fr", "b_to"]
            @assert(haskey(line, key), "line $i: a line without linecode, should specify $key.")
        end
    end
end


function create_line(id, f_bus, t_bus, n_conductors; kwargs...)
    line = Dict{String,Any}()
    line["id"] = id
    line["f_bus"] = f_bus
    line["t_bus"] = t_bus
    line["n_conductors"] = n_conductors
    add_kwarg!(line, kwargs, :f_terminals, collect(1:n_conductors))
    add_kwarg!(line, kwargs, :t_terminals, collect(1:n_conductors))
    add_kwarg!(line, kwargs, :rs, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(line, kwargs, :xs, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(line, kwargs, :g_fr, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(line, kwargs, :b_fr, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(line, kwargs, :g_to, fill(0.0, n_conductors, n_conductors))
    add_kwarg!(line, kwargs, :b_to, fill(0.0, n_conductors, n_conductors))
    return line
end


function create_line_with_linecode(id, f_bus, t_bus, linecode, length; kwargs...)
    line = Dict{String,Any}()
    line["id"] = id
    line["f_bus"] = f_bus
    line["t_bus"] = t_bus
    line["linecode"] = linecode["id"]
    line["length"] = length
    n_conductors = linecode["n_conductors"]
    add_kwarg!(line, kwargs, :f_terminals, collect(1:n_conductors))
    add_kwarg!(line, kwargs, :t_terminals, collect(1:n_conductors))
    return line
end


# Voltage zone

DTYPES[:voltage_zone] = Dict(
    :id => AbstractString,
    :buses => Array{String, 1},
    :vnom => Real,
    :vm_ln_min => Array{<:Real, 1},
    :vm_ln_max => Array{<:Real, 1},
    :vm_lg_min => Array{<:Real, 1},
    :vm_lg_max => Array{<:Real, 1},
    :vm_ng_min => Array{<:Real, 1},
    :vm_ng_max => Array{<:Real, 1},
    :vm_ll_min => Array{<:Real, 1},
    :vm_ll_max => Array{<:Real, 1},
)


CHECKS[:voltage_zone] = function check_voltage_zone(data, bus)
    i = bus["id"]
    @assert(haskey(bus, "vnom"), "Voltage zone $i: a voltage zone should specify a vnom.")
end


function create_voltage_zone(id, vnom; kwargs...)
    voltage_zone = Dict{String, Any}(
        "id" => id,
        "vnom" => vnom,
        "buses" => Array{String, 1}(),
    )
end

# Bus

DTYPES[:bus] = Dict(
    :id => AbstractString,
    :bus => String,
    :bus_type => Int,
    :terminals => Array{<:Int},
    :phases => Array{<:Int},
    :neutral => Union{Int, Missing},
    :grounded => Bool,
    :voltage_zone => String,
    :rg => Real,
    :xg => Real,
    :vnom => Real,
    :vm_ln_min => Array{<:Real, 1},
    :vm_ln_max => Array{<:Real, 1},
    :vm_lg_min => Array{<:Real, 1},
    :vm_lg_max => Array{<:Real, 1},
    :vm_ng_min => Array{<:Real, 1},
    :vm_ng_max => Array{<:Real, 1},
    :vm_ll_min => Array{<:Real, 1},
    :vm_ll_max => Array{<:Real, 1},
    :vm_fix => Array{<:Real, 1},
    :va_fix => Array{<:Real, 1},
)


CHECKS[:bus] = function check_bus(data, bus)
    i = bus["id"]
    if !haskey(bus, "neutral")
        @assert(!haskey(bus, "vm_ln_max"), "Bus $i does not have a neutral, and therefore vm_ln_max bounds do not make sense.")
        @assert(!haskey(bus, "vm_ln_min"), "Bus $i does not have a neutral, and therefore vm_ln_min bounds do not make sense.")
        @assert(!haskey(bus, "vm_ng_min"), "Bus $i does not have a neutral, and therefore a vm_ng_min bound does not make sense.")
        @assert(!haskey(bus, "vm_ng_max"), "Bus $i does not have a neutral, and therefore a vm_ng_max bound does not make sense.")
    end
    nph = length(bus["phases"])
    for key in ["vm_ln_max", "vm_ln_min", "vm_lg_min", "vm_lg_max"]
        if haskey(bus, key)
            @assert(length(bus[key])==nph, "For bus $i, the length of $key should match the number of phases.")
        end
    end
    if nph==1
        @assert(!haskey(bus, "vm_ll_min"), "Bus $i only has a single phase, so line-to-line bounds do not make sense.")
        @assert(!haskey(bus, "vm_ll_max"), "Bus $i only has a single phase, so line-to-line bounds do not make sense.")
    else
        els = nph==2 ? 1 : nph
        for key in ["vm_ll_min", "vm_ll_max"]
            if haskey(bus, key)
                @assert(length(bus["vm_ll_min"])==els, "For bus $i, the $key bound should have only $els elements.")
            end
        end
    end
    if haskey(bus, "voltage_zone")
        zone_id = bus["voltage_zone"]
        @assert(haskey(data["voltage_zone"], zone_id), "Bus $i: the voltage zone $zone_id is not defined.")
        voltage_zone = data["voltage_zone"][zone_id]
        for key in ["vm_ln_min", "vm_ln_max", "vm_lg_min", "vm_lg_max", "vm_ng_min", "vm_ng_max", "vm_ll_min", "vm_ll_max"]
            if haskey(voltage_zone, key)
                @assert(!haskey(bus, key), "Bus $i: the property $key is already specified in the voltage_zone; this cannot be overwritten.")
            end
        end
    end
end

function create_bus(id; kwargs...)
    bus = Dict{String,Any}()
    bus["id"] = id
    add_kwarg!(bus, kwargs, :terminals, collect(1:4))
    add_kwarg!(bus, kwargs, :phases, collect(1:3))
    add_kwarg!(bus, kwargs, :neutral, 4)
    add_kwarg!(bus, kwargs, :grounded, false)
    copy_kwargs_to_dict_if_present!(bus, kwargs, [:vm_nom_ln, :vm_ln_min, :vm_ln_max, :vm_ng_min, :vm_ng_max, :vm_lg_min, :vm_lg_max, :vm_ll_min, :vm_ll_max, :vm_fix, :va_fix, :rg, :xg])
    return bus
end

function create_bus_in_zone(id, voltage_zone, data_model; kwargs...)
    bus = Dict{String,Any}()
    bus["id"] = id
    bus["voltage_zone"] = voltage_zone
    if !(id in data_model["voltage_zone"][voltage_zone]["buses"])
        push!(data_model["voltage_zone"][voltage_zone]["buses"], id)
    end
    add_kwarg!(bus, kwargs, :terminals, collect(1:4))
    add_kwarg!(bus, kwargs, :phases, collect(1:3))
    add_kwarg!(bus, kwargs, :neutral, 4)
    add_kwarg!(bus, kwargs, :grounded, false)
    copy_kwargs_to_dict_if_present!(bus, kwargs, [:vm_nom_ln, :vm_ln_min, :vm_ln_max, :vm_ng_min, :vm_ng_max, :vm_lg_min, :vm_lg_max, :vm_ll_min, :vm_ll_max, :vm_fix, :va_fix, :rg, :xg])
    return bus
end


# Load

DTYPES[:load] = Dict(
    :id => AbstractString,
    :bus => String,
    :terminals => Array{<:Int},
    :configuration => String,
    :pd => Array{<:Real, 1},
    :qd => Array{<:Real, 1},
    :model => String,
    :vnom => Array{<:Real, 1},
    :alpha => Array{<:Real, 1},
    :beta => Array{<:Real, 1},
)


CHECKS[:load] = function check_load(data, load)
    N = length(load["pd"])
    i = load["id"]

    @assert(load["configuration"] in ["wye", "delta"])

    if load["configuration"]=="delta"
        @assert(N>=3,  "Load $i: delta-connected loads should have at least dimension 3, not $N.")
    end

    for key in [:qd, :vm_nom, :alpha, :beta]
        if haskey(load, key)
            @assert(length(load[key])==N, "Load $i: $key should have the same length as pd.")
        end
    end
end


function create_load(id, bus; kwargs...)
    load = Dict{String,Any}()
    load["id"] = id
    load["bus"] = bus

    add_kwarg!(load, kwargs, :configuration, "wye")
    add_kwarg!(load, kwargs, :terminals, load["configuration"]=="wye" ? [1, 2, 3, 4] : [1, 2, 3])
    add_kwarg!(load, kwargs, :pd, fill(0.0, 3))
    add_kwarg!(load, kwargs, :qd, fill(0.0, 3))
    add_kwarg!(load, kwargs, :model, "constant_power")
    copy_kwargs_to_dict_if_present!(load, kwargs, [:vm_nom, :alpha, :beta])
    return load
end

# generatorerator


DTYPES[:generator] = Dict(
    :id => AbstractString,
    :bus => String,
    :terminals => Array{<:Int},
    :configuration => String,
    :pd_set => Array{<:Real, 1},
    :qd_set => Array{<:Real, 1},
    :pd_min => Array{<:Real, 1},
    :pd_max => Array{<:Real, 1},
    :pd_min => Array{<:Real, 1},
    :pd_max => Array{<:Real, 1},
)


CHECKS[:generator] = function check_generator(data, generator)
    N = length(generator["pd_set"])
    i = generator["id"]

    @assert(generator["configuration"] in ["wye", "delta"])

    if load["configuration"]=="delta"
        @assert(N>=3,  "generatorerator $i: delta-connected loads should have at least dimension 3, not $N.")
    end
end


function create_generator(id, bus; kwargs...)
    generator = Dict{String,Any}()
    generator["id"] = id
    generator["bus"] = bus
    add_kwarg!(generator, kwargs, :configuration, "wye")
    add_kwarg!(generator, kwargs, :terminals, generator["configuration"]=="wye" ? [1, 2, 3, 4] : [1, 2, 3])
    copy_kwargs_to_dict_if_present!(generator, kwargs, [:pd_min, :pd_max, :qd_min, :qd_max])
    return generator
end


# Transformer, n-windings three-phase lossy


DTYPES[:transformer_nw_lossy] = Dict(
    :id => AbstractString,
    :n_windings=>Int,
    :n_phases=>Int,
    :buses => Array{Int, 1},
    :terminals => Array{Array{Int, 1}, 1},
    :vnom => Array{<:Real, 1},
    :snom => Array{<:Real, 1},
    :configuration => Array{String, 1},
    :polarity => Array{Bool, 1},
    :xsc => Array{<:Real, 1},
    :rs => Array{<:Real, 1},
    :loadloss => Real,
    :imag => Real,
    :tm_fix => Array{Array{Bool, 1}, 1},
    :tm_set => Array{<:Array{<:Real, 1}, 1},
    :tm_min => Array{<:Array{<:Real, 1}, 1},
    :tm_max => Array{<:Array{<:Real, 1}, 1},
    :tm_step => Array{<:Array{<:Real, 1}, 1},
)


CHECKS[:transformer_nw_lossy] = function check_transformer_nw_lossy(data, trans)
    @assert(all([x in ["wye", "delta"] for x in trans["configuration"]]))
    n_windings = trans["n_windings"]
    @assert(length(trans["xsc"])==n_windings^2-n_windings)
end


function create_transformer_nw_lossy(id, n_windings, buses, terminals, vnom, snom; kwargs...)
    trans = Dict{String,Any}()
    trans["id"] = id
    trans["buses"] = buses
    trans["n_windings"] = n_windings
    trans["terminals"] = terminals
    trans["vnom"] = vnom
    trans["snom"] = snom
    add_kwarg!(trans, kwargs, :configuration, fill("wye", n_windings))
    add_kwarg!(trans, kwargs, :polarity, fill(true, n_windings))
    add_kwarg!(trans, kwargs, :rs, fill(0.0, n_windings))
    add_kwarg!(trans, kwargs, :xsc, fill(0.0, n_windings^2-n_windings))
    add_kwarg!(trans, kwargs, :loadloss, 0.0)
    add_kwarg!(trans, kwargs, :imag, 0.0)
    add_kwarg!(trans, kwargs, :tm_set, fill(fill(1.0, 3), n_windings))
    add_kwarg!(trans, kwargs, :tm_min, fill(fill(0.9, 3), n_windings))
    add_kwarg!(trans, kwargs, :tm_max, fill(fill(1.1, 3), n_windings))
    add_kwarg!(trans, kwargs, :tm_step, fill(fill(1/32, 3), n_windings))
    add_kwarg!(trans, kwargs, :tm_fix, fill(fill(true, 3), n_windings))
    copy_kwargs_to_dict_if_present!(trans, kwargs, [:tm_min, :tm_max])
    return trans
end


# Transformer, two-winding three-phase

DTYPES[:transformer_2w_ideal] = Dict(
    :id => AbstractString,
    :f_bus => String,
    :t_bus => String,
    :configuration => String,
    :f_terminals => Array{Int, 1},
    :t_terminals => Array{Int, 1},
    :tm_nom => Real,
    :tm_set => Real,
    :tm_min => Real,
    :tm_max => Real,
    :tm_step => Real,
    :tm_fix => Real,
)


CHECKS[:transformer_2w_ideal] = function check_transformer_2w_ideal(data, trans)
end


function create_transformer_2w_ideal(id, f_bus, t_bus, tm_nom; kwargs...)
    trans = Dict{String,Any}()
    trans["id"] = id
    trans["f_bus"] = f_bus
    trans["t_bus"] = t_bus
    trans["tm_nom"] = tm_nom
    add_kwarg!(trans, kwargs, :configuration, "wye")
    add_kwarg!(trans, kwargs, :f_terminals, trans["configuration"]=="wye" ? [1, 2, 3, 4] : [1, 2, 3])
    add_kwarg!(trans, kwargs, :t_terminals, [1, 2, 3, 4])
    add_kwarg!(trans, kwargs, :tm_set, 1.0)
    add_kwarg!(trans, kwargs, :tm_min, 0.9)
    add_kwarg!(trans, kwargs, :tm_max, 1.1)
    add_kwarg!(trans, kwargs, :tm_step, 1/32)
    add_kwarg!(trans, kwargs, :tm_fix, true)
    return trans
end


# Capacitor

DTYPES[:capacitor] = Dict(
    :id => AbstractString,
    :bus => String,
    :terminals => Array{Int, 1},
    :configuration => String,
    :qd_ref => Array{<:Real, 1},
    :vm_nom => Real,
)


CHECKS[:capacitor] = function check_capacitor(data, cap)
    i = cap["id"]
    N = length(cap["terminals"])
    config = cap["configuration"]
    if config=="wye"
        @assert(length(cap["qd_ref"])==N-1, "Capacitor $i: qd_ref should have $(N-1) elements.")
    else
        @assert(length(cap["qd_ref"])==N, "Capacitor $i: qd_ref should have $N elements.")
    end
    @assert(config in ["delta", "wye", "wye-grounded", "wye-floating"])
    if config=="delta"
        @assert(N>=3, "Capacitor $i: delta-connected capacitors should have at least 3 elements.")
    end
end


function create_capacitor(id, bus, vm_nom; kwargs...)
    cap = Dict{String,Any}()
    cap["id"] = id
    cap["bus"] = bus
    cap["vm_nom"] = vm_nom
    add_kwarg!(cap, kwargs, :configuration, "wye")
    add_kwarg!(cap, kwargs, :terminals, collect(1:4))
    add_kwarg!(cap, kwargs, :qd_ref, fill(0.0, 3))
    return cap
end


# Shunt

DTYPES[:shunt] = Dict(
    :id => AbstractString,
    :bus => String,
    :terminals => Array{Int, 1},
    :g_sh => Array{<:Real, 2},
    :b_sh => Array{<:Real, 2},
)


CHECKS[:shunt] = function check_shunt(data, shunt)
end


function create_shunt(id, bus, terminals; kwargs...)
    shunt = Dict{String,Any}()
    shunt["id"] = id
    shunt["bus"] = bus
    shunt["terminals"] = terminals
    add_kwarg!(shunt, kwargs, :g_sh, fill(0.0, length(terminals), length(terminals)))
    add_kwarg!(shunt, kwargs, :b_sh, fill(0.0, length(terminals), length(terminals)))
    return shunt
end
