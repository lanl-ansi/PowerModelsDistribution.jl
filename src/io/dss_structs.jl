# Defines data structures (defaults) for OpenDSS objects


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

    Ys = (complex(0.0, 2 * pi * basefreq * c1) * 2.0 + complex(0.0, 2 * pi * basefreq * c0)) / 3.0
    Ym = (complex(0.0, 2 * pi * basefreq * c0) - complex(0.0, 2 * pi * basefreq * c1)) / 3.0

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

        # TODO: convert to rmatrix, xmatrix?
    elseif haskey(kwargs, "rmatrix") && haskey(kwargs, "xmatrix")
        rmatrix = parse_matrix(Float64, get(kwargs, "rmatrix"))
        xmatrix = parse_matrix(Float64, get(kwargs, "xmatrix"))
    elseif haskey(kwargs, "z1")
        z1 = complex(parse_array(Float64, get(kwargs, "z1"))...)
        z2 = complex(parse_array(Float64, get(kwargs, "z2", z1))...)
        z0 = complex(parse_array(Float64, get(kwargs, "z0", z1))...)

        Z = zeros(Complex64, phases, phases)

        for i in 1:phases
            if phases == 1
                Z[i,i] = complex(z1...) / 3.0
            else
                Z[i,i] = (complex(z2...) + complex(z1...) + complex(z0...)) / 3.0
            end
        end

        if phases == 3
            Z[2, 1] = Z[3, 2] = Z[1, 3] = (conj(exp(-2*pi*im/3))^2 * z2 + conj(exp(-2*pi*im/3)) * z1 + z0) / 3
            Z[3, 1] = Z[1, 2] = Z[2, 3] = (conj(exp(-2*pi*im/3))^2 * z1 + conj(exp(-2*pi*im/3)) * z2 + z0) / 3
        end

        rmatrix = real(Z)
        xmatrix = imag(Z)
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
                               "z1" => [real(z1), imag(z1)],
                               "z2" => [real(z2), imag(z2)],
                               "z0" => [real(z0) imag(z0)],
                               "z" => [real(z), imag(z)],
                               "rcurve" => get(kwargs, "rcurve", ""),
                               "lcurve" => get(kwargs, "lcurve", ""),
                               "lmh" => lmh,
                               # Inherited Properties
                               "normamps" => get(kwargs, "normamps", 400.0),
                               "emergamps" => get(kwargs, "emergamps", 600.0),
                               "repair" => get(kwargs, "repair", 3.0),
                               "faultrate" => get(kwargs, "faultrate", 0.1),
                               "pctperm" => get(kwargs, "pctperm", 20.0),
                               "basefreq" => get(kwargs, "basefreq", 60.0),
                               "enabled" => get(kwargs, "enabled", false)
                              )
    return reactor
end


""
function createVSource(bus1::Int, name::String, bus2::Int=0; kwargs...)
    x1r1 = get(kwargs, "x1r1", 4.0)
    x0r0 = get(kwargs, "x0r0", 3.0)

    basekv = get(kwargs, "basekv", 115.0)
    pu = get(kwargs, "pu", 1.0)
    rs = 0.0
    rm = 0.0
    xs = 0.1
    xm = 0.0

    phases = get(kwargs, "phases", 3)
    factor = phases == 1 ? 1.0 : sqrt(3.0)

    mvasc3 = get(kwargs, "mvasc3", 2000.0)
    mvasc1 = get(kwargs, "mvasc1", 2100.0)

    isc3 = get(kwargs, "isc3", 10000.0)
    isc1 = get(kwargs, "isc1", 10500.0)

    r1 = get(kwargs, "r1", 1.65)
    x1 = get(kwargs, "x1", 6.6)
    r0 = get(kwargs, "r0", 1.9)
    x0 = get(kwargs, "x0", 5.7)
    r2 = r1
    x2 = x1

    z1 = parse_array(Float64, get(kwargs, "z1", [0.0, 0.0]))
    z2 = parse_array(Float64, get(kwargs, "z2", [0.0, 0.0]))
    z0 = parse_array(Float64, get(kwargs, "z0", [0.0, 0.0]))

    puz1 = parse_array(Float64, get(kwargs, "puz1", [0.0, 0.0]))
    puz2 = parse_array(Float64, get(kwargs, "puz2", [0.0, 0.0]))
    puz0 = parse_array(Float64, get(kwargs, "puz0", [0.0, 0.0]))

    basemva = get(kwargs, "basemva", 100.0)

    Zbase = basekv^2 / basemva

    if (haskey(kwargs, "mvasc3") && haskey(kwargs, "mvasc1")) || (haskey(kwargs, "isc3") && haskey(isc1, "isc1"))
        if haskey(kwargs, "mvasc3") && haskey(kwargs, "mvasc1")
            mvasc3 = get(kwargs, "mvasc3")
            mvasc1 = get(kwargs, "mvasc1")

            isc3 = mvasc3 * 1e3 / (basekv * sqrt(3.0))
            isc1 = mvasc1 * 1e3 / (basekv * factor)
        elseif haskey(kwargs, "isc3") && haskey(isc1, "isc1")
            isc3 = get(kwargs, "isc3")
            isc1 = get(kwargs, "isc1")

            mvasc3 = sqrt(3) * basekv * isc3 / 1e3
            mvasc1 = factor * basekv * isc1 / 1e3
        end

        x1 = basekv^2 / mvasc3 / sqrt(1.0 + 1.0 / x1r1^2)

        r1 = x1 / x1r
        r2 = r1
        x2 = x1

        a = 1.0 + x0r0^2
        b = 4.0*(r1 + x1 * x0r0)
        c = 4.0 * (r1^2 + x1^2)- (3.0 * basekv * 1000.0 / factor / Isc1)^2
        r0 = max((-b + sqrt(b^2 - 4 * a * c)) / (2 * a), (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
        x0 = r0 * x0r0

        xs = (2.0 * x1 + x0) / 3.0
        rs = (2.0 * r1 + r0) / 3.0

        rm = (r0 - r1) / 3.0
        xm = (x0 - x1) / 3.0
    elseif any([haskey(kwargs, key) for key in ["r1", "x1", "z1", "puz1"]])
        if haskey(kwargs, "puz1")
            puz1 = complex(parse_array(Float64, get(kwargs, "puz1"))...)
            puz2 = complex(parse_array(Float64, get(kwargs, "puz2", puz1))...)
            puz0 = complex(parse_array(Float64, get(kwargs, "puz0", puz1))...)

            r1 = real(puz1) * Zbase
            x1 = imag(puz1) * Zbase
            r2 = real(puz2) * Zbase
            x2 = imag(puz2) * Zbase
            r0 = real(puz0) * Zbase
            x1 = imag(puz0) * Zbase
        elseif (haskey(kwargs, "r1") && haskey(kwargs, "x1"))
            r1 = get(kwargs, "r1")
            x1 = get(kwargs, "x1")

            r2 = get(kwargs, "r2", r1)
            x2 = get(kwargs, "x2", x1)
            r0 = get(kwargs, "r0", r1)
            x0 = get(kwargs, "x0", x1)
        elseif haskey(kwargs, "z1")
            z1 = complex(parse_array(Float64, get(kwargs, "z1")))
            z2 = complex(parse_array(Float64, get(kwargs, "z2"), z1))
            z0 = complex(parse_array(Float64, get(kwargs, "z0"), z1))

            r1 = real(z1)
            x1 = imag(z1)
            r2 = real(z2)
            x2 = imag(z2)
            r0 = real(z0)
            x0 = imag(z0)
        end

        isc3 = basekv * 1e3 / sqrt(3.0) * abs(complex(r1, x1))

        if phases == 1
            r0 = r1
            x0 = x1
            r2 = r1
            x2 = x1
        end

        rs = (2.0 * r1 + r0) / 3.0
        xs = (2.0 * x1 + x0) / 3.0

        isc1 = basekv * 1e3 / factor / abs(complex(rs, xs))

        mvasc3 = sqrt(3) * basekv * isc3 / 1e3
        mvasc1 = factor * basekv * isc1 / 1e3

        xm = xs - x1

        rs = (2.0 * r1 + r0) / 3.0
        rm = (r0 - r1) / 3.0
    end

    Z = zeros(Complex64, phases, phases)
    if r1 == r2 && x1 == x2
        Zs = complex(Rs, Xs)
        Zm = complex(Rm, Xm)

        for i in 1:phases
            Z[i,i] = Zs
            for j in 1:i-1
                Z[i, j] = Z[j, i] = Zm
            end
        end
    else
        z1 = complex(r1, x1)
        z2 = complex(r2, x2)
        z0 = complex(r0, x0)

        for i in 1:phases
            Z[i,i] = (z1 + z2 + z0) / 3.0
        end

        if phases == 3
            Z[2, 1] = Z[3, 2] = Z[1, 3] = (conj(exp(-2*pi*im/3))^2 * z2 + conj(exp(-2*pi*im/3)) * z1 + z0) / 3
            Z[3, 1] = Z[1, 2] = Z[2, 3] = (conj(exp(-2*pi*im/3))^2 * z1 + conj(exp(-2*pi*im/3)) * z2 + z0) / 3
        end
    end

    Vmag = phases == 1 ? basekv * pu * 1e3 : basekv * pu * 1e3 / 2 / sin(pi / phases)

    if !haskey(kwargs, "puz1") && Zbase > 0.0
        puz1 = complex(r1 / Zbase, x1 / Zbase)
        puz2 = complex(r2 / Zbase, x2 / Zbase)
        puz0 = complex(r0 / Zbase, x0 / Zbase)
    end

    vsource = Dict{String,Any}("bus1" => bus1,
                               "basekv" => basekv,
                               "pu" => pu,
                               "angle" => get(kwargs, "angle", 0.0),
                               "frequency" => get(kwargs, "frequency"),
                               "phases" => phases,
                               "mvasc3" => mvasc3,
                               "mvasc1" => mvasc1,
                               "x1r1" => x1r1,
                               "x0r0" => x0r0,
                               "isc3" => isc3,
                               "isc1" => isc1,
                               "r1" => r1,
                               "x1" => x1,
                               "r0" => r0,
                               "x0" => x0,
                               "scantype" => get(kwargs, "scantype", "pos"),
                               "sequence" => get(kwargs, "sequence", "pos"),
                               "bus2" => bus2,
                               "z1" => [real(z1), imag(z1)],
                               "z0" => [real(z0), imag(z0)],
                               "z2" => [real(z2), imag(z2)],
                               "puz1" => [real(puz1), imag(puz1)],
                               "puz0" => [real(puz0), imag(puz0)],
                               "puz2" => [real(puz2), imag(puz2)],
                               "basemva" => basemva,
                               "yearly" => parse_array(Float64, get(kwargs, "yearly", get(kwargs, "daily", [1.0, 1.0]))),
                               "daily" => parse_array(Float64, get(kwargs, "daily", [1.0, 1.0])),
                               "duty" => get(kwargs, "duty", ""),
                               # Inherited Properties
                               "spectrum" => get(kwargs, "spectrum", "defaultvsource"),
                               "basefreq" => get(kwargs, "basefreq", 60.0),
                               "enabled" => get(kwargs, "enabled", false)
                               # Derived Properties
                               "rmatrix" => real(Z),
                               "xmatrix" => imag(Z),
                               "vmag" => Vmag,
                              )

    return vsource
end


""
function createTransformer(bus1, bus2, name, bus3::Int=0; kwargs...)
    
    transformer = Dict{String,Any}()

    return transformer
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
