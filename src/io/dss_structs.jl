# Defines data structures (defaults) for OpenDSS objects


""
function createLoad(bus1::Int, name::AbstractString; kwargs...)
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
                            )

    return load
end


""
function createLine(bus1, bus2, name; kwargs...)
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

    rmatrix = parse_matrix(get(kwargs, "rmatrix", real(Z)))
    xmatrix = parse_matrix(get(kwargs, "xmatrix", imag(Z)))
    cmatrix = parse_matrix(get(kwargs, "cmatrix", imag(Yc)))

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
                            "enabled" => true
                            )

    return line
end
