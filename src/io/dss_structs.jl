# Defines data structures (defaults) for OpenDSS objects


""
function parse_conn(conn::String)::String
    if conn in ["wye", "y", "ln"]
        return "wye"
    elseif conn in ["delta", "ll"]
        return "delta"
    else
        warn(LOGGER, "Unsupported connection $conn, defaulting to \"wye\"")
        return "wye"
    end
end


""
function createLoad(bus1::Int, name::String; kwargs...)
    kv = get(kwargs, "kv", 12.47)
    kw = get(kwargs, "kw", 10.0)
    pf = get(kwargs, "pf", 0.88)
    kvar = get(kwargs, "kvar", 5.0)
    kva = get(kwargs, "kva", kw / pf)

    if haskey(kwargs, "kw") && haskey(kwargs, "pf")
        kvar = sign(pf) * kw * sqrt(1.0 / pf^2 - 1.0)
        kva = abs(kw) + kvar^2
    elseif haskey(kwargs, "kw") && haskey(kwargs, "kvar")
        kva = abs(kw) + kvar^2
        if kva > 0.0
            pf = kw / kva
            if kvar != 0.0
                pf *= sign(kw * kvar)
            end
        end
    elseif haskey(kwargs, "kva") && haskey(kwargs, "pf")
        kw = kva * abs(pf)
        kvar = sign(pf) * kw * sqrt(1.0 / pf^2 - 1.0)
    elseif haskey(kwargs, "pf") && kwargs["pf"] != 0.88
            kvar = sign(pf) * kw * sqrt(1.0 / pf^2 - 1.0)
            kva = abs(kw) + kvar^2
    end

    # TODO: yearly, daily, duty, growth, model
    # TODO: ZIPV (7 coefficient array, depends on model keyword)

    if haskey(kwargs, "like")
        warn(LOGGER, "\"like\" keyword on load $name is not supported.")
    end

    load = Dict{String,Any}("name" => name,
                            "phases" => get(kwargs, "phases", 3),
                            "bus1" => bus1,
                            "kv" => kv,
                            "kw" => kw,
                            "pf" => pf,
                            "model" => get(kwargs, "model", 1),
                            "yearly" => get(kwargs, "yearly", get(kwargs, "daily", complex(1.0, 1.0))),
                            "daily" => get(kwargs, "daily", complex(1.0, 1.0)),
                            "duty" => get(kwargs, "duty", ""),
                            "growth" => get(kwargs, "growth", ""),
                            "conn" => get(kwargs, "conn", "wye"),
                            "kvar" => kvar,
                            "rneut" => get(kwargs, "rneut", -1.0),
                            "xneut" => get(kwargs, "xneut", 0.0),
                            "status" => get(kwargs, "status", "variable"),
                            "class" => get(kwargs, "class", 1),
                            "vminpu" => get(kwargs, "vminpu", 0.95),
                            "vmaxpu" => get(kwargs, "vmaxpu", 1.05),
                            "vminnorm" => get(kwargs, "vminnorm", 0.0),
                            "vminemerg" => get(kwargs, "vminemerg", 0.0),
                            "xfkva" => get(kwargs, "xfkva", 0.0),
                            "allocationfactor" => get(kwargs, "allocationfactor", 0.5),
                            "kva" => kva
                            "%mean" => get(kwargs, "%mean", 0.5),
                            "%stddev" => get(kwargs, "%stddev", 0.1),
                            "cvrwatts" => get(kwargs, "cvrwatts", 1.0),
                            "cvrvars" => get(kwargs, "cvrvars", 2.0),
                            "kwh" => get(kwargs, "kwh", 0.0),
                            "kwhdays" => get(kwargs, "kwhdays", 30.0),
                            "cfactor" => get(kwargs, "cfactor", 4.0),
                            "cvrcurve" => get(kwargs, "cvrcurve", ""),
                            "numcust" => get(kwargs, "numcust", 1),
                            "zipv" => get(kwargs, "zipv", ""),
                            "%seriesrl" => get(kwargs, "%seriesrl", 0.5),
                            "relweight" => get(kwargs, "relweight", 1.0),
                            "vlowpu" => get(kwargs, "vlowpu", 0.5),
                            "puxharm" => get(kwargs, "puxharm", 0.0),
                            "xrharm" => get(kwargs, "xrharm", 6.0),
                            # Inherited Properties
                            "spectrum" => get(kwargs, "spectrum", "defaultload"),
                            "basefreq" => get(kwargs, "basefreq", 60.0)
                            "enabled" => get(kwargs, "enabled", true)
                           )

    return load
end


""
function createLine(bus1::Int, bus2::Int, name::String; kwargs...)
    phases = get(kwargs, "phases", 3)
    basefreq = get(kwargs, "basefreq", 60.0)

    r1 = get(kwargs, "r1", 0.058)
    x1 = get(kwargs, "x1", 0.1206)
    c1 = get(kwargs, "c1", 3.4)

    Ztemp = complex(r1, x1) * 2.0

    if phases == 1
        r0 = r1
        x0 = x1
        c0 = c1
    else
        r0 = get(kwargs, "r0", 0.1784)
        x0 = get(kwargs, "x0", 0.4047)
        c0 = get(kwargs, "c0", 1.6)
    end

    c1 = haskey(kwargs, "b1") ? get(kwargs, "b1") / (2 * pi * basefreq) * 1.0e-6 : c1
    c0 = haskey(kwargs, "b0") ? get(kwargs, "b0") / (2 * pi * basefreq) * 1.0e-6 : c0

    Zs = (Ztemp + complex(r0, x0)) / 3.0
    Zm = (complex(r0, x0) - complex(r1, x1)) / 3.0

    Yc1 = 2 * pi * basefreq * c1
    Yc0 = 2 * pi * basefreq * c0

    Ys = (complex(0., Yc1) * 2.0 + complex(0.0, Yc0)) / 3.0
    Ym = (complex(0.0, Yc0) - complex(0.0, Yc1)) / 3.0

    Z = zeros(Complex64, phases, phases)
    Yc = zeros(Complex64, phases, phases)
    for i in 1:phases
        Z[i,i] = Zs
        Yc[i,i] = Ys
        for j in 1:i-1
            Z[i,j] = Z[j,i] = Zm
            Yc[i,j] = Yc[j,i] = Ym
        end
    end

    rmatrix = parse_matrix(Float64, get(kwargs, "rmatrix", real(Z)))
    xmatrix = parse_matrix(Float64, get(kwargs, "xmatrix", imag(Z)))
    cmatrix = parse_matrix(Float64, get(kwargs, "cmatrix", imag(Yc)))

    # TODO: rg, xg
    rg = get(kwargs, "rg", 0.01805)
    xg = get(kwargs, "xg", 0.155081)

    # TODO: switch
    if get(kwargs, "switch", false)
        warn(LOGGER, "\"switch\" keyword in line $name is not supported.")
    end

    if haskey(kwargs, "like")
        warn(LOGGER, "\"like\" keyword on load $name is not supported.")
    end

    line = Dict{String,Any}("name" => name,
                            "bus1" => bus1,
                            "bus2" => bus2,
                            "linecode" => get(kwargs, "linecode", ""),
                            "length" => get(kwargs, "length", 1.0),
                            "phases" => phases,
                            "r1" => r1,
                            "x1" => x1,
                            "r0" => r0,
                            "x0" => x0,
                            "c1" => c1,
                            "c0" => c0,
                            "rmatrix" => rmatrix,
                            "xmatrix" => xmatrix,
                            "cmatrix" => cmatrix,
                            "switch" => get(kwargs, "switch", false)
                            "rg" => rg,
                            "xg" => xg,
                            "rho" => get(kwargs, "rho", 100),
                            "geometry" => get(kwargs, "geometry"),
                            "units" => get(kwargs, "units", "none"),
                            "spacing" => get(kwargs, "spacing", ""),
                            "wires" => get(kwargs, "wires", ""),
                            "earthmodel" => get(kwargs, "earthmodel", ""),
                            "cncables" => get(kwargs, "cncables", ""),
                            "tscables" => get(kwargs, "tscables", ""),
                            "b1" => b1,
                            "b0" => b0,
                            # Inherited Properties
                            "normamps" => get(kwargs, "normamps", 400.0)
                            "emergamps" => get(kwargs, "emergamps", 600.0)
                            "faultrate" => get(kwargs, "faultrate", 0.1)
                            "pctperm" => get(kwargs, "pctperm", 20.0)
                            "repair" => get(kwargs, "repair", 3.0)
                            "basefreq" => basefreq
                            "enabled" => get(kwargs, "enabled", true)
                           )

    return line
end


""
function createGenerator(bus1::Int, name::String; kwargs...)
    conn = parse_conn(get(kwargs, "conn", "wye"))

    kw = get(kwargs, "kw", 100.0)
    kva = get(kwargs, "kva", kw * 1.2)
    kvar = get(kwargs, "kvar", 60.0)
    kvarmax = get(kwargs, "maxkvar", kvar * 2.0)
    kvarmin = get(kwargs, "minkvar", -kvarmax)

    gen = Dict{String,Any}("name" => name,
                           "phases" => get(kwargs, "phases", 3),
                           "bus1" => bus1,
                           "kv" => get(kwargs, "kv", 12.47),
                           "kw" => kw,
                           "pf" => get(kwargs, "pf", 0.80),
                           "model" => get(kwargs, "model", 1),
                           "yearly" => get(kwargs, "yearly", get(kwargs, "daily", complex(0.0, 0.0))),
                           "daily" => get(kwargs, "daily", get(kwargs, "daily", complex(0.0, 0.0))),
                           "duty" => get(kwargs, "duty", ""),
                           "dispmode" => get(kwargs, "dispmode", "default"),
                           "disvalue" => get(kwargs, "dispvalue", 0.0),
                           "conn" => conn,
                           "kvar" => kvar,
                           "rneut" => get(kwargs, "rneut", 0.0),
                           "xneut" => get(kwargs, "xneut", 0.0),
                           "status" => get(kwargs, "status", "variable"),
                           "class" => get(kwargs, "class", 1),
                           "vpu" => get(kwargs, "vpu", 1.0),
                           "maxkvar" => kvarmax,
                           "minkvar" => kvarmin,
                           "pvfactor" => get(kwargs, "pvfactor", 0.1),
                           "debugtrace" => get(kwargs, "debugtrace", false),
                           "vminpu" => get(kwargs, "", 0.9),
                           "vmaxpu" => get(kwargs, "", 1.10),
                           "forceon" => get(kwargs, "forceon", false),
                           "kva" => kva,
                           "mva" => kva * 0.001,
                           "xd" => get(kwargs, "xd", 1.0),
                           "xdp" => get(kwargs, "xdp", 0.28),
                           "xdpp" => get(kwargs, "xdpp", 0.20),
                           "h" => get(kwargs, "h", 1.0),
                           "d" => get(kwargs, "d", 1.0),
                           "usermodel" => get(kwargs, "usermodel", ""),
                           "userdata" => get(kwargs, "userdata", ""),
                           "shaftmodel" => get(kwargs, "shaftmodel", ""),
                           "shaftdata" => get(kwargs, "shaftdata", ""),
                           "dutystart" => get(kwargs, "dutystart", 0.0),
                           "balanced" => get(kwargs, "balanced", false),
                           "xrdp" => get(kwargs, "xrdp", 20.0),
                           # Inherited Properties
                           "spectrum" => get(kwargs, "spectrum", "defaultgen"),
                           "basefreq" => get(kwargs, "basefreq", 60.0),
                           "enabled" => get(kwargs, "enabled", true)
                          )
    return gen
end


""
function createCapacitor(bus1::Int, name::String, bus2::Int=0; kwargs...)
    phases = get(kwargs, "phases", 3)

    capacitor = Dict{String,Any}("bus1" => bus1,
                                 "bus2" => bus2,
                                 "phases" => phases,
                                 "kvar" => get(kwargs, "kvar", 1200.0),
                                 "kv" => get(kwargs, "kv", 12.47),
                                 "conn" => parse_conn(get(kwargs, "conn", "wye")),
                                 "cmatrix" => parse_matrix(Float64, get(kwargs, "cmatrix", zeros(phases, phases))),
                                 "cuf" => parse_array(Float64, get(kwargs, "cuf", zeros(phases))),
                                 "r" => parse_array(Float64, get(kwargs, "r", zeros(phases))),
                                 "xl" => parse_array(Float64, get(kwargs, "xl", zeros(phases))),
                                 "harm" => parse_array(Float64, get(kwargs, "harm", zeros(phases))),
                                 "numsteps" => get(kwargs, "numsteps", 1),
                                 "states" => parse_array(Bool, get(kwargs, "states", zeros(Bool, phases))),
                                 # Inherited Properties
                                 "normamps" => get(kwargs, "normamps", 400.0)
                                 "emergamps" => get(kwargs, "emergamps", 600.0)
                                 "faultrate" => get(kwargs, "faultrate", 0.1)
                                 "pctperm" => get(kwargs, "pctperm", 20.0)
                                 "basefreq" => get(kwargs, "basefreq", 60.0)
                                 "enabled" => get(kwargs, "enabled", true)
                                )
    return capacitor
end


""
function createReactor(bus1::Int, name::String, bus2::Int=0; kwargs...)
    phases = get(kwargs, "phases", 3)
    kvar = get(kwargs, "kvar", 1200.0)
    kv = get(kwargs, "kv", 12.47)
    conn = parse_conn(get(kwargs, "conn", "wye"))
    parallel = get(kwargs, "parallel", false)

    normamps = get(kwargs, "normamps", 400.0)
    emergamps = get(kwargs, "emergamps", 600.0)

    rp = get(kwargs, "rp", 0.0)

    if (haskey(kwargs, "kv") && haskey(kwargs, "kvar")) || haskey(kwargs, "x") || haskey(kwargs, "lmh") || haskey(kwargs, "z")
        r = get(kwargs, "r")

        if haskey(kwargs, "kvar") && haskey("kv")
            kvarperphase = kvar / phases
            if conn == "delta"
                phasekv = kv
            else
                if phases == 2 || phases == 3
                    phasekv = kv / sqrt(3.0)
                else
                    phasekv = kv
                end
            end

            x = phasekv^2 * 1.0e3 / kvarperphase
            l = x / (2 * pi) / basefreq
            normamps = kvarperphase / phasekv
            emergamps = normamps * 1.35

        elseif haskey(kwargs, "x")
            x = get(kwargs, "x")
            l = x / (2 * pi) / basefreq

        elseif haskey(kwargs, "lmh")
            l = get(kwargs, "lmh") / 1.0e3
            x = l * 2 * pi * basefreq

        elseif haskey(kwargs, "z")
            z = complex(parse_array(Float64, get(kwargs, "z")...))
            r = real(z)
            x = imag(z)
            l = x / (2 * pi) / basefreq
        end

        Value = inv(complex(r, l * 2 * pi * FYpriFreq))

        if rp != 0.0
            Caccum(Value, complex(Gp, 0.0))
        end

        Value2 = Value * 2.0

        Value = -Value

        YPrimTemp = zeros(phases, phases)

        if conn == "delta"
            for i in 1:phases
                YPrimTemp[i,i] = Value2
                for j in 1:i-1
                    YPrimTemp[i,j] = YPrimTemp[j,i] = Value
                end
            end
        else
            for i in 1:phases
                YPrimTemp[i,i] = Value
            end
        end




    elseif haskey(kwargs, "rmatrix") && haskey(kwargs, "xmatrix")
        rmatrix = parse_matrix(Float64, get(kwargs, "rmatrix"))
        xmatrix = parse_matrix(Float64, get(kwargs, "xmatrix"))
    elseif haskey(kwargs, "z1")
        z1 = get(kwargs, "z1")
        z2 = get(kwargs, "z2", z1)
        z0 = get(kwargs, "z0", z1)
    else
        warn(LOGGER, "Reactor $name is not adequately defined")
    end

    reactor = Dict{String,Any}("bus1" => bus1,
                               "bus2" => bus2,
                               "phases" => phases,
                               "kvar" => kvar
                               "kv" => kv
                               "conn" = conn,
                               "rmatrix" => rmatrix,
                               "xmatrix" => xmatrix,
                               "parallel" => parallel
                               "r" => r,
                               "x" => x,
                               "rp" => rp,
                               "z1" => z1,
                               "z2" => z2,
                               "z0" => z0,
                               "z" => z,
                               "rcurve" => get(kwargs, "rcurve", ""),
                               "lcurve" => get(kwargs, "lcurve", ""),
                               "lmh" => lmh,
                               # Inherited Properties
                               "normamps" => get(kwargs, "normamps", 400.0),
                               "emergamps" => get(kwargs, "emergamps", 600.0),
                               "repair" = get(kwargs, "repair", 3.0),
                               "faultrate" => get(kwargs, "faultrate", 0.1),
                               "pctperm" => get(kwargs, "pctperm", 20.0),
                               "basefreq" => get(kwargs, "basefreq", 60.0),
                               "enabled" => get(kwargs, "enabled", false)
                              )
    return reactor
end

""
function get_dtypes(comp::String)::Dict
    default_dicts = Dict{String,Any}("line" => createLine(0, 0, ""),
                                     "load" => createLoad(0, "")
                                     "generator" => createGenerator(0, "")
                                     "capacitor" => createCapacitor(0, "")
                                     "reactor" => createReactor(0, "")
                                     "transformer" => createTransformer(0, 0, "")
                                     "linecode" => createLinecode(0, 0, "")
                                     "circuit" => createCircuit("")
                                    )

    return Dict{String,Type}((k, typeof(v) for (k, v) in default_dicts[comp]))
end

get_dtypes(comp::String, key::String)::Type = get_dtypes(comp)[key]
