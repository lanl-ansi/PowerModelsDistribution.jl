# Defines data structures (defaults) for OpenDSS objects
import LinearAlgebra: diagm
import Statistics: mean, std


const _convert_to_meters = Dict{String,Float64}(
    "mi" => 1609.3,
    "km" => 1000.0,
    "kft" => 304.8,
    "m" => 1.0,
    "ft" => 0.3048,
    "in" => 0.0254,
    "cm" => 0.01,
    "mm" => 0.001,
    "none" => 1.0
)


"""
Creates a Dict{String,Any} containing all of the properties of a Linecode. See
OpenDSS documentation for valid fields and ways to specify the different
properties.
"""
function _create_linecode(name::String=""; kwargs...)::Dict{String,Any}
    phases = get(kwargs, :nphases, 3)
    circuit_basefreq = get(kwargs, :circuit_basefreq, 60.0)
    basefreq = get(kwargs, :basefreq, circuit_basefreq)

    r1 = get(kwargs, :r1, 0.058)
    x1 = get(kwargs, :x1, 0.1206)

    if haskey(kwargs, :b1)
        b1 = kwargs[:b1]
        c1 = b1 / (2 * pi * basefreq)
    else
        c1 = get(kwargs, :c1, 3.4)
        b1 = c1 * (2 * pi * basefreq)
    end

    if phases == 1
        r0 = r1
        x0 = x1
        c0 = c1
        b0 = b1
    else
        r0 = get(kwargs, :r0, 0.1784)
        x0 = get(kwargs, :x0, 0.4047)
        if haskey(kwargs, :b0)
            b0 = kwargs[:b0]
            c0 = b0 / (2 * pi * basefreq)
        else
            c0 = get(kwargs, :c0, 1.6)
            b0 = c0 * (2 * pi * basefreq)
        end
    end

    Zs = (complex(r1, x1) * 2.0 + complex(r0, x0)) / 3.0
    Zm = (complex(r0, x0) - complex(r1, x1)) / 3.0

    Ys = (complex(0.0, 2 * pi * basefreq * c1) * 2.0 + complex(0.0, 2 * pi * basefreq * c0)) / 3.0
    Ym = (complex(0.0, 2 * pi * basefreq * c0) - complex(0.0, 2 * pi * basefreq * c1)) / 3.0

    Z = zeros(Complex{Float64}, phases, phases)
    Yc = zeros(Complex{Float64}, phases, phases)
    for i in 1:phases
        Z[i,i] = Zs
        Yc[i,i] = Ys
        for j in 1:i-1
            Z[i,j] = Z[j,i] = Zm
            Yc[i,j] = Yc[j,i] = Ym
        end
    end

    rmatrix = get(kwargs, :rmatrix, real(Z))
    xmatrix = get(kwargs, :xmatrix, imag(Z))
    cmatrix = get(kwargs, :cmatrix, imag(Yc) / (2 * pi * basefreq))

    units = get(kwargs, :units, "none")

    Dict{String,Any}(
        "name" => name,
        "nphases" => phases,
        "r1" => r1 / _convert_to_meters[units],
        "x1" => x1 / _convert_to_meters[units],
        "r0" => r0 / _convert_to_meters[units],
        "x0" => x0 / _convert_to_meters[units],
        "c1" => c1 / _convert_to_meters[units],
        "c0" => c0 / _convert_to_meters[units],
        "units" => "m",
        "rmatrix" => rmatrix / _convert_to_meters[units],
        "xmatrix" => xmatrix / _convert_to_meters[units],
        "cmatrix" => cmatrix / _convert_to_meters[units],
        "basefreq" => basefreq,
        "normamps" => get(kwargs, :normamps, 400.0),
        "emergamps" => get(kwargs, :emergamps, 600.0),
        "faultrate" => get(kwargs, :faultrate, 0.1),
        "pctperm" => get(kwargs, :pctperm, 20.0),
        "repair" => get(kwargs, :repair, 3.0),
        "kron" => get(kwargs, :kron, false),
        "rg" => get(kwargs, :rg, 0.01805),
        "xg" => get(kwargs, :xg, 0.155081),
        "rho" => get(kwargs, :rho, 100.0),
        "neutral" => get(kwargs, :neutral, 3),
        "b1" => b1 / _convert_to_meters[units],
        "b0" => b0 / _convert_to_meters[units],
        "like" => get(kwargs, :like, "")
    )
end


"""
Creates a Dict{String,Any} containing all of the properties for a Line. See
OpenDSS documentation for valid fields and ways to specify the different
properties.
"""
function _create_line(name::String=""; kwargs...)::Dict{String,Any}
    bus1 = get(kwargs, :bus1, "")
    bus2 = get(kwargs, :bus2, "")

    phases = get(kwargs, :phases, 3)
    circuit_basefreq = get(kwargs, :circuit_basefreq, 60.0)
    basefreq = get(kwargs, :basefreq, circuit_basefreq)

    r1 = get(kwargs, :r1, 0.058)
    x1 = get(kwargs, :x1, 0.1206)

    if haskey(kwargs, :b1)
        b1 = kwargs[:b1]
        c1 = b1 / (2 * pi * basefreq)
    else
        c1 = get(kwargs, :c1, 3.4)
        b1 = c1 * (2 * pi * basefreq)
    end

    if phases == 1
        r0 = r1
        x0 = x1
        c0 = c1
        b0 = b1
    else
        r0 = get(kwargs, :r0, 0.1784)
        x0 = get(kwargs, :x0, 0.4047)
        if haskey(kwargs, :b0)
            b0 = kwargs[:b0]
            c0 = b0 / (2 * pi * basefreq)
        else
            c0 = get(kwargs, :c0, 1.6)
            b0 = c0 * (2 * pi * basefreq)
        end
    end

    Zs = (complex(r1, x1) * 2.0 + complex(r0, x0)) / 3.0
    Zm = (complex(r0, x0) - complex(r1, x1)) / 3.0

    Ys = (complex(0.0, 2 * pi * basefreq * c1) * 2.0 + complex(0.0, 2 * pi * basefreq * c0)) / 3.0
    Ym = (complex(0.0, 2 * pi * basefreq * c0) - complex(0.0, 2 * pi * basefreq * c1)) / 3.0

    Z  = Matrix{Complex{Float64}}(undef, phases, phases)
    Yc = Matrix{Complex{Float64}}(undef, phases, phases)
    for i in 1:phases
        Z[i,i] = Zs
        Yc[i,i] = Ys
        for j in 1:i-1
            Z[i,j] = Z[j,i] = Zm
            Yc[i,j] = Yc[j,i] = Ym
        end
    end

    rmatrix = get(kwargs, :rmatrix, real(Z))
    xmatrix = get(kwargs, :xmatrix, imag(Z))
    cmatrix = get(kwargs, :cmatrix, imag(Yc) / (2 * pi * basefreq))

    rg = get(kwargs, :rg, 0.01805)
    xg = get(kwargs, :xg, 0.155081)
    rho = get(kwargs, :rho, 100.0)

    # TODO: support length mismatch between line and linecode?
    # Currently this does not change the values of rmatrix and xmatrix due to
    # lenmult=1, code is only in place for future use.
    lenmult = 1.0

    kxg = xg / log(658.5 * sqrt(rho / circuit_basefreq))
    xgmod = xg != 0.0 ?  0.5 * kxg * log(basefreq / circuit_basefreq) : 0.0

    units = get(kwargs, :units, "none")
    len = get(kwargs, :switch, false) ? 0.001 : get(kwargs, :length, 1.0) * _convert_to_meters[units]

    if haskey(kwargs, :rg)
        Memento.warn(_LOGGER, "Rg,Xg are not fully supported")
    end

    rmatrix .+= rg * (basefreq / circuit_basefreq - 1.0)
    rmatrix .*= lenmult
    xmatrix .-= xgmod
    xmatrix .*= lenmult * (basefreq / circuit_basefreq)

    Dict{String,Any}(
        "name" => name,
        "bus1" => bus1,
        "bus2" => bus2,
        "linecode" => get(kwargs, :linecode, ""),
        "length" => len,
        "phases" => phases,
        "r1" => r1 / _convert_to_meters[units],
        "x1" => x1 / _convert_to_meters[units],
        "r0" => r0 / _convert_to_meters[units],
        "x0" => x0 / _convert_to_meters[units],
        "c1" => c1 / _convert_to_meters[units],
        "c0" => c0 / _convert_to_meters[units],
        "rmatrix" => rmatrix / _convert_to_meters[units],
        "xmatrix" => xmatrix / _convert_to_meters[units],
        "cmatrix" => cmatrix / _convert_to_meters[units],
        "switch" => get(kwargs, :switch, false),
        "rg" => rg / _convert_to_meters[units],
        "xg" => xg / _convert_to_meters[units],
        "rho" => get(kwargs, :rho, 100),
        "geometry" => get(kwargs, :geometry, ""),
        "units" => "m",
        "spacing" => get(kwargs, :spacing, ""),
        "wires" => get(kwargs, :wires, Vector{String}([])),
        "earthmodel" => get(kwargs, :earthmodel, ""),
        "cncables" => get(kwargs, :cncables, Vector{String}([])),
        "tscables" => get(kwargs, :tscables, Vector{String}([])),
        "b1" => b1 / _convert_to_meters[units],
        "b0" => b0 / _convert_to_meters[units],
        # Inherited Properties
        "normamps" => get(kwargs, :normamps, 400.0),
        "emergamps" => get(kwargs, :emergamps, 600.0),
        "faultrate" => get(kwargs, :faultrate, 0.1),
        "pctperm" => get(kwargs, :pctperm, 20.0),
        "repair" => get(kwargs, :repair, 3.0),
        "basefreq" => basefreq,
        "enabled" => get(kwargs, :enabled, true),
        "like" => get(kwargs, :like, "")
    )
end


"""
Creates a Dict{String,Any} containing all of the expected properties for a
Load. See OpenDSS documentation for valid fields and ways to specify the
different properties.
"""
function _create_load(name::String=""; kwargs...)::Dict{String,Any}
    bus1 = get(kwargs, :bus1, "")

    kv = get(kwargs, :kv, 12.47)
    kw = get(kwargs, :kw, 10.0)
    pf = get(kwargs, :pf, 0.88)
    kvar = get(kwargs, :kvar, 5.0)
    kva = get(kwargs, :kva, kw / pf)

    if haskey(kwargs, :kw) && haskey(kwargs, :pf)
        kvar = sign(pf) * kw * sqrt(1.0 / pf^2 - 1.0)
        kva = abs(kw) + kvar^2
    elseif haskey(kwargs, :kw) && haskey(kwargs, :kvar)
        kva = abs(kw) + kvar^2
        if kva > 0.0
            pf = kw / kva
            if kvar != 0.0
                pf *= sign(kw * kvar)
            end
        end
    elseif haskey(kwargs, :kva) && haskey(kwargs, :pf)
        kw = kva * abs(pf)
        kvar = sign(pf) * kw * sqrt(1.0 / pf^2 - 1.0)
    elseif haskey(kwargs, :pf) && kwargs[:pf] != 0.88
            kvar = sign(pf) * kw * sqrt(1.0 / pf^2 - 1.0)
            kva = abs(kw) + kvar^2
    end

    # TODO: ZIPV (7 coefficient array, depends on model keyword)

    Dict{String,Any}(
        "name" => name,
        "phases" => get(kwargs, :phases, 3),
        "bus1" => bus1,
        "kv" => kv,
        "kw" => kw,
        "pf" => pf,
        "model" => get(kwargs, :model, 1),
        "yearly" => get(kwargs, :yearly, ""),
        "daily" => get(kwargs, :daily, ""),
        "duty" => get(kwargs, :duty, ""),
        "growth" => get(kwargs, :growth, ""),
        "conn" => get(kwargs, :conn, WYE),
        "kvar" => kvar,
        "rneut" => get(kwargs, :rneut, -1.0),
        "xneut" => get(kwargs, :xneut, 0.0),
        "status" => get(kwargs, :status, "variable"),
        "class" => get(kwargs, :class, 1),
        "vminpu" => get(kwargs, :vminpu, 0.95),
        "vmaxpu" => get(kwargs, :vmaxpu, 1.05),
        "vminnorm" => get(kwargs, :vminnorm, 0.0),
        "vminemerg" => get(kwargs, :vminemerg, 0.0),
        "xfkva" => get(kwargs, :xfkva, 0.0),
        "allocationfactor" => get(kwargs, :allocationfactor, 0.5),
        "kva" => kva,
        "%mean" => get(kwargs, Symbol("%mean"), 0.5),
        "%stddev" => get(kwargs, Symbol("%stddev"), 0.1),
        "cvrwatts" => get(kwargs, :cvrwatts, 1.0),
        "cvrvars" => get(kwargs, :cvrvars, 2.0),
        "kwh" => get(kwargs, :kwh, 0.0),
        "kwhdays" => get(kwargs, :kwhdays, 30.0),
        "cfactor" => get(kwargs, :cfactor, 4.0),
        "cvrcurve" => get(kwargs, :cvrcurve, ""),
        "numcust" => get(kwargs, :numcust, 1),
        "zipv" => get(kwargs, :zipv, ""),
        "%seriesrl" => get(kwargs, Symbol("%seriesrl"), 0.5),
        "relweight" => get(kwargs, :relweight, 1.0),
        "vlowpu" => get(kwargs, :vlowpu, 0.5),
        "puxharm" => get(kwargs, :puxharm, 0.0),
        "xrharm" => get(kwargs, :xrharm, 6.0),
        # Inherited Properties
        "spectrum" => get(kwargs, :spectrum, "defaultload"),
        "basefreq" => get(kwargs, :basefreq, 60.0),
        "enabled" => get(kwargs, :enabled, true),
        "like" => get(kwargs, :like, "")
    )
end


"""
Creates a Dict{String,Any} containing all of the expected properties for a
Generator. See OpenDSS documentation for valid fields and ways to specify the
different properties.
"""
function _create_generator(name::String=""; kwargs...)::Dict{String,Any}
    bus1 = get(kwargs, :bus1, "")

    conn = get(kwargs, :conn, WYE)

    kw = get(kwargs, :kw, 100.0)
    kva = get(kwargs, :kva, kw * 1.2)
    kvar = get(kwargs, :kvar, 60.0)
    kvarmax = get(kwargs, :maxkvar, kvar * 2.0)
    kvarmin = get(kwargs, :minkvar, -kvarmax)

    Dict{String,Any}(
        "name" => name,
        "phases" => get(kwargs, :phases, 3),
        "bus1" => bus1,
        "kv" => get(kwargs, :kv, 12.47),
        "kw" => kw,
        "pf" => get(kwargs, :pf, 0.80),
        "model" => get(kwargs, :model, 1),
        "yearly" => get(kwargs, :yearly, get(kwargs, :daily, Vector{Float64}([0.0, 0.0]))),
        "daily" => get(kwargs, :daily, Vector{Float64}([0.0, 0.0])),
        "duty" => get(kwargs, :duty, ""),
        "dispmode" => get(kwargs, :dispmode, "default"),
        "disvalue" => get(kwargs, :dispvalue, 0.0),
        "conn" => conn,
        "kvar" => kvar,
        "rneut" => get(kwargs, :rneut, 0.0),
        "xneut" => get(kwargs, :xneut, 0.0),
        "status" => get(kwargs, :status, "variable"),
        "class" => get(kwargs, :class, 1),
        "vpu" => get(kwargs, :vpu, 1.0),
        "maxkvar" => kvarmax,
        "minkvar" => kvarmin,
        "pvfactor" => get(kwargs, :pvfactor, 0.1),
        "debugtrace" => get(kwargs, :debugtrace, false),
        "vminpu" => get(kwargs, :vminpu, 0.9),
        "vmaxpu" => get(kwargs, :vmaxpu, 1.10),
        "forceon" => get(kwargs, :forceon, false),
        "kva" => kva,
        "mva" => kva * 0.001,
        "xd" => get(kwargs, :xd, 1.0),
        "xdp" => get(kwargs, :xdp, 0.28),
        "xdpp" => get(kwargs, :xdpp, 0.20),
        "h" => get(kwargs, :h, 1.0),
        "d" => get(kwargs, :d, 1.0),
        "usermodel" => get(kwargs, :usermodel, ""),
        "userdata" => get(kwargs, :userdata, ""),
        "shaftmodel" => get(kwargs, :shaftmodel, ""),
        "shaftdata" => get(kwargs, :shaftdata, ""),
        "dutystart" => get(kwargs, :dutystart, 0.0),
        "balanced" => get(kwargs, :balanced, false),
        "xrdp" => get(kwargs, :xrdp, 20.0),
        # Inherited Properties
        "spectrum" => get(kwargs, :spectrum, "defaultgen"),
        "basefreq" => get(kwargs, :basefreq, 60.0),
        "enabled" => get(kwargs, :enabled, true),
        "like" => get(kwargs, :like, "")
    )
end


"""
Creates a Dict{String,Any} containing all of the expected properties for a
Capacitor. If `bus2` is not specified, the capacitor will be treated as a shunt.
See OpenDSS documentation for valid fields and ways to specify the
different properties.
"""
function _create_capacitor(name::String=""; kwargs...)::Dict{String,Any}
    bus1 = get(kwargs, :bus1, "")

    phases = get(kwargs, :phases, 3)

    bus2 = get(kwargs, :bus2, string(split(bus1, ".")[1],".",join(fill("0", phases), ".")))

    Dict{String,Any}(
        "name" => name,
        "bus1" => bus1,
        "bus2" => bus2,
        "phases" => phases,
        "kvar" => get(kwargs, :kvar, 1200.0),
        "kv" => get(kwargs, :kv, 12.47),
        "conn" => get(kwargs, :conn, WYE),
        "cmatrix" => get(kwargs, :cmatrix, zeros(phases, phases)),
        "cuf" => get(kwargs, :cuf, zeros(phases)),
        "r" => get(kwargs, :r, zeros(phases)),
        "xl" => get(kwargs, :xl, zeros(phases)),
        "harm" => get(kwargs, :harm, zeros(phases)),
        "numsteps" => get(kwargs, :numsteps, 1),
        "states" => get(kwargs, :states, zeros(Bool, phases)),
        # Inherited Properties
        "normamps" => get(kwargs, :normamps, 400.0),
        "emergamps" => get(kwargs, :emergamps, 600.0),
        "faultrate" => get(kwargs, :faultrate, 0.1),
        "pctperm" => get(kwargs, :pctperm, 20.0),
        "basefreq" => get(kwargs, :basefreq, 60.0),
        "enabled" => get(kwargs, :enabled, true),
        "like" => get(kwargs, :like, "")
    )
end


"""
Creates a Dict{String,Any} containing all of the expected properties for a
Reactor. If `bus2` is not specified Reactor is treated like a shunt. See
OpenDSS documentation for valid fields and ways to specify the different
properties.
"""
function _create_reactor(name::String=""; kwargs...)::Dict{String,Any}
    bus1 = get(kwargs, :bus1, "")
    bus2 = get(kwargs, :bus2, "")

    phases = get(kwargs, :phases, 3)
    kvar = get(kwargs, :kvar, 1200.0)
    kv = get(kwargs, :kv, 12.47)
    conn = get(kwargs, :conn, WYE)
    parallel = get(kwargs, :parallel, false)

    normamps = get(kwargs, :normamps, 400.0)
    emergamps = get(kwargs, :emergamps, 600.0)
    basefreq = get(kwargs, :basefreq, 60.0)

    r = get(kwargs, :r, 0.0)
    x = get(kwargs, :x, abs(kv*1e3))

    rp = get(kwargs, :rp, 0.0)

    z1 = get(kwargs, :z1, [r, x])
    z0 = get(kwargs, :z0, [0.0, 0.0])
    z2 = get(kwargs, :z2, [0.0, 0.0])
    z = get(kwargs, :z, [r, x])

    lmh = get(kwargs, :lmh, x / (2 * pi * basefreq) * 1e3)

    # TODO: handle `parallel`
    if (haskey(kwargs, :kv) && haskey(kwargs, :kvar)) || haskey(kwargs, :x) || haskey(kwargs, :lmh) || haskey(kwargs, :z)
        if haskey(kwargs, :kvar) && haskey(kwargs, :kv)
            kvarperphase = kvar / phases
            if conn == DELTA
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
            lmh = l * 1e3
            normamps = kvarperphase / phasekv
            emergamps = normamps * 1.35

        elseif haskey(kwargs, :x)
            x = kwargs[:x]
            l = x / (2 * pi) / basefreq
            lmh = l * 1e3

        elseif haskey(kwargs, :lmh)
            l = kwargs[:lmh] / 1.0e3
            x = l * 2 * pi * basefreq

        elseif haskey(kwargs, :z)
            z = complex(kwargs[:z]...)
            r = real(z)
            x = imag(z)
            l = x / (2 * pi) / basefreq
            lmh = l * 1e3
        end

        rmatrix = diagm(0 => fill(r, phases))
        xmatrix = diagm(0 => fill(x, phases))
    elseif haskey(kwargs, :rmatrix) && haskey(kwargs, :xmatrix)
        rmatrix = kwargs[:rmatrix]
        xmatrix = kwargs[:xmatrix]

        r = rmatrix[1, 1]
        x = xmatrix[1, 1]

        # TODO: account for off-diagonal and single phase
        z = z1 = z2 = z0 = complex(r, x)
        lmh = x / (2 * pi) / basefreq * 1e3
    elseif haskey(kwargs, :z1)
        z1 = complex(kwargs[:z1]...)
        z2 = complex(get(kwargs, :z2, z1)...)
        z0 = complex(get(kwargs, :z0, z1)...)

        Z = zeros(Complex{Float64}, phases, phases)

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

        r = rmatrix[1,1]
        x = xmatrix[1,1]
        lmh = x / (2 * pi * basefreq) * 1e3
    else
        rmatrix = diagm(0 => fill(r, phases))
        xmatrix = diagm(0 => fill(x, phases))
    end

    Dict{String,Any}(
        "name" => name,
        "bus1" => bus1,
        "bus2" => bus2,
        "phases" => phases,
        "kvar" => kvar,
        "kv" => kv,
        "conn" => conn,
        "rmatrix" => rmatrix,
        "xmatrix" => xmatrix,
        "parallel" => parallel,
        "r" => r,
        "x" => x,
        "rp" => rp,
        "z1" => [real(z1), imag(z1)],
        "z2" => [real(z2), imag(z2)],
        "z0" => [real(z0) imag(z0)],
        "z" => [real(z), imag(z)],
        "rcurve" => get(kwargs, :rcurve, ""),
        "lcurve" => get(kwargs, :lcurve, ""),
        "lmh" => lmh,
        # Inherited Properties
        "normamps" => normamps,
        "emergamps" => emergamps,
        "repair" => get(kwargs, :repair, 3.0),
        "faultrate" => get(kwargs, :faultrate, 0.1),
        "pctperm" => get(kwargs, :pctperm, 20.0),
        "basefreq" => basefreq,
        "enabled" => get(kwargs, :enabled, true),
        "like" => get(kwargs, :like, "")
    )
end


"""
Creates a Dict{String,Any} containing all of the expected properties for a
Voltage Source. If `bus2` is not specified, VSource will be treated like a
generator. Mostly used as `source` which represents the circuit. See
OpenDSS documentation for valid fields and ways to specify the different
properties.
"""
function _create_vsource(name::String=""; kwargs...)::Dict{String,Any}
    phases = get(kwargs, :phases, 3)

    bus1 = get(kwargs, :bus1, "sourcebus.$(join(1:phases, "."))")
    bus2 = get(kwargs, :bus2, replace(bus1, r"\.\d" => ".0"))

    x1r1 = get(kwargs, :x1r1, 4.0)
    x0r0 = get(kwargs, :x0r0, 3.0)

    basekv = get(kwargs, :basekv, 115.0)
    pu = get(kwargs, :pu, 1.0)
    rs = 0.0
    rm = 0.0
    xs = 0.1
    xm = 0.0

    factor = phases == 1 ? 1.0 : sqrt(3.0)

    mvasc3 = get(kwargs, :mvasc3, 2000.0)
    mvasc1 = get(kwargs, :mvasc1, 2100.0)

    isc3 = get(kwargs, :isc3, 10041.0)
    isc1 = get(kwargs, :isc1, 10543.0)

    r1 = get(kwargs, :r1, 1.65)
    x1 = get(kwargs, :x1, 6.6)
    r0 = get(kwargs, :r0, 1.9)
    x0 = get(kwargs, :x0, 5.7)
    r2 = r1
    x2 = x1

    z1 = get(kwargs, :z1, [0.0, 0.0])
    z2 = get(kwargs, :z2, [0.0, 0.0])
    z0 = get(kwargs, :z0, [0.0, 0.0])

    puz1 = get(kwargs, :puz1, [0.0, 0.0])
    puz2 = get(kwargs, :puz2, [0.0, 0.0])
    puz0 = get(kwargs, :puz0, [0.0, 0.0])

    basemva = get(kwargs, :basemva, 100.0)

    Zbase = basekv^2 / basemva

    if (haskey(kwargs, :mvasc3) || haskey(kwargs, :mvasc1)) || (haskey(kwargs, :isc3) || haskey(kwargs, :isc1))
        if haskey(kwargs, :mvasc3) || haskey(kwargs, :mvasc1)
            isc3 = haskey(kwargs, :mvasc3) ? mvasc3 * 1e3 / (basekv * sqrt(3.0)) : isc3
            isc1 = haskey(kwargs, :mvasc1) ? mvasc1 * 1e3 / (basekv * factor) : isc1
        elseif haskey(kwargs, :isc3) || haskey(kwargs, :isc1)
            mvasc3 = haskey(kwargs, :isc3) ? sqrt(3) * basekv * isc3 / 1e3 : mvasc3
            mvasc1 = haskey(kwargs, :isc1) ? factor * basekv * isc1 / 1e3 : mvasc1
        end

        x1 = basekv^2 / mvasc3 / sqrt(1.0 + 1.0 / x1r1^2)

        r1 = x1 / x1r1
        r2 = r1
        x2 = x1

        a = 1.0 + x0r0^2
        b = 4.0*(r1 + x1 * x0r0)
        c = 4.0 * (r1^2 + x1^2)- (3.0 * basekv * 1000.0 / factor / isc1)^2
        r0 = max((-b + sqrt(b^2 - 4 * a * c)) / (2 * a), (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
        x0 = r0 * x0r0

        xs = (2.0 * x1 + x0) / 3.0
        rs = (2.0 * r1 + r0) / 3.0

        rm = (r0 - r1) / 3.0
        xm = (x0 - x1) / 3.0
    elseif any([haskey(kwargs, key) for key in [:r1, :x1, :z1, :puz1]])
        if haskey(kwargs, :puz1)
            puz1 = complex(kwargs[:puz1]...)
            puz2 = complex(get(kwargs, :puz2, puz1)...)
            puz0 = complex(get(kwargs, :puz0, puz1)...)

            r1 = real(puz1) * Zbase
            x1 = imag(puz1) * Zbase
            r2 = real(puz2) * Zbase
            x2 = imag(puz2) * Zbase
            r0 = real(puz0) * Zbase
            x1 = imag(puz0) * Zbase
        elseif (haskey(kwargs, :r1) && haskey(kwargs, :x1))
            r1 = kwargs[:r1]
            x1 = kwargs[:x1]

            r2 = get(kwargs, :r2, r1)
            x2 = get(kwargs, :x2, x1)
            r0 = get(kwargs, :r0, r1)
            x0 = get(kwargs, :x0, x1)
        elseif haskey(kwargs, :z1)
            z1 = complex(kwargs[:z1]...)
            z2 = complex(get(kwargs, :z2, z1)...)
            z0 = complex(get(kwargs, :z0, z1)...)

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

        rm = (r0 - r1) / 3.0
        xm = xs - x1

        isc1 = basekv * 1e3 / factor / abs(complex(rs, xs))

        mvasc3 = sqrt(3) * basekv * isc3 / 1e3
        mvasc1 = factor * basekv * isc1 / 1e3
    else
        rs = (2.0 * r1 + r0) / 3.0
        xs = (2.0 * x1 + x0) / 3.0

        rm = (r0 - r1) / 3.0
        xm = (x0 - x1) / 3.0
    end

    Z = zeros(Complex{Float64}, phases, phases)
    if r1 == r2 && x1 == x2
        Zs = complex(rs, xs)
        Zm = complex(rm, xm)

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

    if !haskey(kwargs, :puz1) && Zbase > 0.0
        puz1 = complex(r1 / Zbase, x1 / Zbase)
        puz2 = complex(r2 / Zbase, x2 / Zbase)
        puz0 = complex(r0 / Zbase, x0 / Zbase)
    end

    Dict{String,Any}(
        "name" => name,
        "bus1" => bus1,
        "basekv" => basekv,
        "pu" => pu,
        "angle" => get(kwargs, :angle, 0.0),
        "frequency" => get(kwargs, :frequency, get(kwargs, :basefreq, 60.0)),
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
        "scantype" => get(kwargs, :scantype, "pos"),
        "sequence" => get(kwargs, :sequence, "pos"),
        "bus2" => bus2,
        "z1" => [real(z1), imag(z1)],
        "z0" => [real(z0), imag(z0)],
        "z2" => [real(z2), imag(z2)],
        "puz1" => [real(puz1), imag(puz1)],
        "puz0" => [real(puz0), imag(puz0)],
        "puz2" => [real(puz2), imag(puz2)],
        "basemva" => basemva,
        "yearly" => get(kwargs, :yearly, get(kwargs, :daily, [1.0, 1.0])),
        "daily" => get(kwargs, :daily, [1.0, 1.0]),
        "duty" => get(kwargs, :duty, ""),
        # Inherited Properties
        "spectrum" => get(kwargs, :spectrum, "defaultvsource"),
        "basefreq" => get(kwargs, :basefreq, 60.0),
        "enabled" => get(kwargs, :enabled, true),
        # Derived Properties
        "rmatrix" => real(Z),
        "xmatrix" => imag(Z),
        "vmag" => Vmag,
        "like" => get(kwargs, :like, "")
    )
end


"alias _create_circuit to _create_vsource"
_create_circuit = _create_vsource


"""
Creates a Dict{String,Any} containing all of the expected properties for a
Transformer. See OpenDSS documentation for valid fields and ways to specify the
different properties.
"""
function _create_transformer(name::String=""; kwargs...)
    windings = isempty(name) ? 3 : get(kwargs, :windings, 2)
    phases = get(kwargs, :phases, 3)

    prcnt_rs = fill(0.2, windings)
    if haskey(kwargs, Symbol("%rs"))
        prcnt_rs = kwargs[Symbol("%rs")]
    elseif haskey(kwargs, Symbol("%loadloss"))
        prcnt_rs[1] = prcnt_rs[2] = kwargs[Symbol("%loadloss")] / 2.0
    end

    temp = Dict{String,Any}("buss" => get(kwargs, :buses, fill("", windings)),
                            "taps" => get(kwargs, :taps, fill(1.0, windings)),
                            "conns" => get(kwargs, :conns, fill(WYE, windings)),
                            "kvs" => get(kwargs, :kvs, fill(12.47, windings)),
                            "kvas" => get(kwargs, :kvas, fill(10.0, windings)),
                            "%rs" => prcnt_rs,
                            "rneuts" => fill(0.0, windings),
                            "xneuts" => fill(0.0, windings)
                           )

    for wdg in [:wdg, :wdg_2, :wdg_3]
        if haskey(kwargs, wdg)
            suffix = kwargs[wdg] == 1 ? "" : "_$(kwargs[wdg])"
            for key in [:bus, :tap, :conn, :kv, :kva, Symbol("%r"), :rneut, :xneut]
                subkey = Symbol(string(key, suffix))
                if haskey(kwargs, subkey)
                    temp[string(key, "s")][kwargs[wdg]] = kwargs[subkey]
                end
            end
        end
    end

    trfm = Dict{String,Any}(
        "name" => name,
        "phases" => phases,
        "windings" => windings,
        # Per wdg
        "wdg" => 1,
        "bus" => temp["buss"][1],
        "conn" => temp["conns"][1],
        "kv" => temp["kvs"][1],
        "kva" => temp["kvas"][1],
        "tap" => temp["taps"][1],
        "%r" => temp["%rs"][1],
        "rneut" => temp["rneuts"][1],
        "xneut" => temp["xneuts"][1],

        "wdg_2" => 2,
        "bus_2" => temp["buss"][2],
        "conn_2" => temp["conns"][2],
        "kv_2" => temp["kvs"][2],
        "kva_2" => temp["kvas"][2],
        "tap_2" => temp["taps"][2],
        "%r_2" => temp["%rs"][2],
        "rneut_2" => temp["rneuts"][2],
        "xneut_2" => temp["xneuts"][2],

        # General
        "buses" => temp["buss"],
        "conns" => temp["conns"],
        "kvs" => temp["kvs"],
        "kvas" => temp["kvas"],
        "taps" => temp["taps"],
        "xhl" => get(kwargs, :xhl, 7.0),
        "xht" => get(kwargs, :xht, 35.0),
        "xlt" => get(kwargs, :xlt, 30.0),
        "xscarray" => get(kwargs, :xscarry, ""),
        "thermal" => get(kwargs, :thermal, 2.0),
        "n" => get(kwargs, :n, 0.8),
        "m" => get(kwargs, :m, 0.8),
        "flrise" => get(kwargs, :flrise, 65.0),
        "hsrise" => get(kwargs, :hsrise, 15.0),
        "%loadloss" => get(kwargs, Symbol("%loadloss"), sum(temp["%rs"][1:2])),
        "%noloadloss" => get(kwargs, Symbol("%noloadloss"), 0.0),
        "normhkva" => get(kwargs, :normhkva, 1.1 * temp["kvas"][1]),
        "emerghkva" => get(kwargs, :emerghkva, 1.5 * temp["kvas"][1]),
        "sub" => get(kwargs, :sub, false),
        "maxtap" => get(kwargs, :maxtap, 1.10),
        "mintap" => get(kwargs, :mintap, 0.90),
        "numtaps" => get(kwargs, :numtaps, 32),
        "subname" => get(kwargs, :subname, ""),
        "%imag" => get(kwargs, Symbol("%imag"), 0.0),
        "ppm_antifloat" => get(kwargs, :ppm_antifloat, 1.0),
        "%rs" => temp["%rs"],
        "bank" => get(kwargs, :bank, ""),
        "xfmrcode" => get(kwargs, :xfmrcode, ""),
        "xrconst" => get(kwargs, :xrconst, false),
        "x12" => get(kwargs, :xhl, 7.0),
        "x13" => get(kwargs, :xht, 35.0),
        "x23" => get(kwargs, :xlt, 30.0),
        "leadlag" => get(kwargs, :leadlag, "lag"),
        # Inherited Properties
        "faultrate" => get(kwargs, :faultrate, 0.1),
        "basefreq" => get(kwargs, :basefreq, 60.0),
        "enabled" => get(kwargs, :enabled, true),
        "like" => get(kwargs, :like, "")
    )

    if windings == 3
        trfm3 = Dict{String,Any}(
            "wdg_3" => 3,
            "bus_3" => temp["buss"][3],
            "conn_3" => temp["conns"][3],
            "kv_3" => temp["kvs"][3],
            "kva_3" => temp["kvas"][3],
            "tap_3" => temp["taps"][3],
            "%r_3" => temp["%rs"][3],
            "rneut_3" => temp["rneuts"][3],
            "xneut_3" => temp["xneuts"][3],
        )

        merge!(trfm, trfm3)
    end

    return trfm
end


"Transformer codes contain all of the same properties as a transformer except bus, buses, bank, xfmrcode"
function _create_xfmrcode(name::String=""; kwargs...)
    windings = isempty(name) ? 3 : get(kwargs, :windings, 2)
    phases = get(kwargs, :phases, 3)

    prcnt_rs = fill(0.2, windings)
    if haskey(kwargs, Symbol("%rs"))
        prcnt_rs = kwargs[Symbol("%rs")]
    elseif haskey(kwargs, Symbol("%loadloss"))
        prcnt_rs[1] = prcnt_rs[2] = kwargs[Symbol("%loadloss")] / 2.0
    end

    temp = Dict{String,Any}(
        "taps" => get(kwargs, :taps, fill(1.0, windings)),
        "conns" => get(kwargs, :conns, fill(WYE, windings)),
        "kvs" => get(kwargs, :kvs, fill(12.47, windings)),
        "kvas" => get(kwargs, :kvas, fill(10.0, windings)),
        "%rs" => prcnt_rs,
        "rneuts" => fill(0.0, windings),
        "xneuts" => fill(0.0, windings)
    )

    for wdg in [:wdg, :wdg_2, :wdg_3]
        if haskey(kwargs, wdg)
            suffix = kwargs[wdg] == 1 ? "" : "_$(kwargs[wdg])"
            for key in [:bus, :tap, :conn, :kv, :kva, Symbol("%r"), :rneut, :xneut]
                subkey = Symbol(string(key, suffix))
                if haskey(kwargs, subkey)
                    temp[string(key, "s")][kwargs[wdg]] = kwargs[subkey]
                end
            end
        end
    end

    xfmrcode = Dict{String,Any}(
        "name" => name,
        "phases" => phases,
        "windings" => windings,
        # Per wdg
        "wdg" => 1,
        "conn" => temp["conns"][1],
        "kv" => temp["kvs"][1],
        "kva" => temp["kvas"][1],
        "tap" => temp["taps"][1],
        "%r" => temp["%rs"][1],
        "rneut" => temp["rneuts"][1],
        "xneut" => temp["xneuts"][1],

        "wdg_2" => 2,
        "conn_2" => temp["conns"][2],
        "kv_2" => temp["kvs"][2],
        "kva_2" => temp["kvas"][2],
        "tap_2" => temp["taps"][2],
        "%r_2" => temp["%rs"][2],
        "rneut_2" => temp["rneuts"][2],
        "xneut_2" => temp["xneuts"][2],

        # General
        "conns" => temp["conns"],
        "kvs" => temp["kvs"],
        "kvas" => temp["kvas"],
        "taps" => temp["taps"],
        "xhl" => get(kwargs, :xhl, 7.0),
        "xht" => get(kwargs, :xht, 35.0),
        "xlt" => get(kwargs, :xlt, 30.0),
        "xscarray" => get(kwargs, :xscarry, ""),
        "thermal" => get(kwargs, :thermal, 2.0),
        "n" => get(kwargs, :n, 0.8),
        "m" => get(kwargs, :m, 0.8),
        "flrise" => get(kwargs, :flrise, 65.0),
        "hsrise" => get(kwargs, :hsrise, 15.0),
        "%loadloss" => get(kwargs, Symbol("%loadloss"), sum(temp["%rs"][1:2])),
        "%noloadloss" => get(kwargs, Symbol("%noloadloss"), 0.0),
        "normhkva" => get(kwargs, :normhkva, 1.1 * temp["kvas"][1]),
        "emerghkva" => get(kwargs, :emerghkva, 1.5 * temp["kvas"][1]),
        "sub" => get(kwargs, :sub, false),
        "maxtap" => get(kwargs, :maxtap, 1.10),
        "mintap" => get(kwargs, :mintap, 0.90),
        "numtaps" => get(kwargs, :numtaps, 32),
        "subname" => get(kwargs, :subname, ""),
        "%imag" => get(kwargs, Symbol("%imag"), 0.0),
        "ppm_antifloat" => get(kwargs, :ppm_antifloat, 1.0),
        "%rs" => temp["%rs"],
        "xrconst" => get(kwargs, :xrconst, false),
        "x12" => get(kwargs, :xhl, 7.0),
        "x13" => get(kwargs, :xht, 35.0),
        "x23" => get(kwargs, :xlt, 30.0),
        "leadlag" => get(kwargs, :leadlag, "lag"),
        # Inherited Properties
        "faultrate" => get(kwargs, :faultrate, 0.1),
        "basefreq" => get(kwargs, :basefreq, 60.0),
        "enabled" => get(kwargs, :enabled, true),
        "like" => get(kwargs, :like, "")
    )

    if windings == 3
        xfmrcode3 = Dict{String,Any}(
            "wdg_3" => 3,
            "conn_3" => temp["conns"][3],
            "kv_3" => temp["kvs"][3],
            "kva_3" => temp["kvas"][3],
            "tap_3" => temp["taps"][3],
            "%r_3" => temp["%rs"][3],
            "rneut_3" => temp["rneuts"][3],
            "xneut_3" => temp["xneuts"][3],
        )

        merge!(xfmrcode, xfmrcode3)
    end

    return xfmrcode
end


"""
Creates a Dict{String,Any} containing all of the expected properties for a
PVSystem. See OpenDSS document
https://github.com/tshort/OpenDSS/blob/master/Doc/OpenDSS%20PVSystem%20Model.doc
for valid fields and ways to specify the different properties.
"""
function _create_pvsystem(name::String=""; kwargs...)
    bus1 = get(kwargs, :bus1, "")

    kv = get(kwargs, :kv, 12.47)
    kw = get(kwargs, :kw, 10.0)
    pf = get(kwargs, :pf, 0.88)
    kvar = get(kwargs, :kvar, 5.0)
    kva = get(kwargs, :kva, kw / pf)


    if haskey(kwargs, :kw) && haskey(kwargs, :pf)
        kvar = sign(pf) * kw * sqrt(1.0 / pf^2 - 1.0)
        kva = abs(kw) + kvar^2
    elseif haskey(kwargs, :kw) && haskey(kwargs, :kvar)
        kva = abs(kw) + kvar^2
        if kva > 0.0
            pf = kw / kva
            if kvar != 0.0
                pf *= sign(kw * kvar)
            end
        end
    elseif haskey(kwargs, :kva) && haskey(kwargs, :pf)
        kw = kva * abs(pf)
        kvar = sign(pf) * kw * sqrt(1.0 / pf^2 - 1.0)
    elseif haskey(kwargs, :pf) && kwargs[:pf] != 0.88
            kvar = sign(pf) * kw * sqrt(1.0 / pf^2 - 1.0)
            kva = abs(kw) + kvar^2
    end

    if haskey(kwargs, :like)
        Memento.warn(_LOGGER, "\"like\" keyword on pvsystem $name is not supported.")
    end

    Dict{String,Any}(
        "name" => name,
        "phases" => get(kwargs, :phases, 3),
        "bus1" => bus1,
        "kv" => kv,
        "kw" => kw,
        "pf" => pf,
        "model" => get(kwargs, :model, 1),
        "yearly" => get(kwargs, :yearly, get(kwargs, :daily, [1.0, 1.0])),
        "daily" => get(kwargs, :daily, [1.0, 1.0]),
        "duty" => get(kwargs, :duty, ""),
        "irradiance" => get(kwargs, :irradiance, 0),
        "pmpp" => get(kwargs, :pmpp, 0),
        "temperature" => get(kwargs, :temperature, 0),
        "conn" => get(kwargs, :conn, WYE),
        "kvar" => kvar,
        "kva" => kva,
        "%cutin" => get(kwargs, :cutin, 0), #TODO not sure what to do with this
        "%cutout" => get(kwargs, :cutout, 0), #TODO not sure what to do with this
        "effcurve" => get(kwargs, :effcurve, ""),
        "p-tcurve" => get(kwargs, :ptcurve, ""),
        "%r" => get(kwargs, :r, 0),
        "%x" => get(kwargs, :x, 0.50),
        "vminpu" => get(kwargs, :vminpu, 0.9),
        "vmaxpu" => get(kwargs, :vmaxpu, 1.1),
        "tyearly" => get(kwargs, :tyearly, 0),
        "tduty" => get(kwargs, :tduty, 0),
        "class" => get(kwargs, :class, 0),
        "usermodel" => get(kwargs, :usermodel, ""),
        "userdata" => get(kwargs, :userdata, ""),
        "debugtrace" => get(kwargs, :debugtrace, "no"),
        "spectrum" => get(kwargs, :spectrum, "defaultpvsystem"),
        # Inherited Properties
        "basefreq" => get(kwargs, :basefreq, 60.0),
        "enabled" => get(kwargs, :enabled, true),
        "like" => get(kwargs, :like, "")
    )
end


"""
Creates a Dict{String,Any} containing all expected properties for a storage
element. See OpenDSS documentation for valid fields and ways to specify the
different properties.
"""
function _create_storage(name::String=""; kwargs...)
    Dict{String,Any}(
        "name" => name,
        "%charge" => get(kwargs, :charge, 100.0),
        "%discharge" => get(kwargs, :discharge, 100.0),
        "%effcharge" => get(kwargs, :effcharge, 90.0),
        "%effdischarge" => get(kwargs, :effdischarge, 90.0),
        "%idlingkvar" => get(kwargs, :idlingkvar, 0.0),
        "%idlingkw" => get(kwargs, :idlingkw, 1.0),
        "%r" => get(kwargs, :r, 0.0),
        "%reserve" => get(kwargs, :reserve, 20.0),
        "%stored" => get(kwargs, :stored, 100.0),
        "%x" => get(kwargs, :x, 50.0),
        "basefreq" => get(kwargs, :basefreq, 60.0),
        "bus1" => get(kwargs, :bus1, ""),
        "chargetrigger" => get(kwargs, :chargetrigger, 0.0),
        "class" => get(kwargs, :class, 0),
        "conn" => get(kwargs, :conn, WYE),
        "daily" => get(kwargs, :daily, [1.0, 1.0]),
        "debugtrace" => get(kwargs, :debugtrace, false),
        "dischargetrigger" => get(kwargs, :dischargetrigger, 0.0),
        "dispmode" => get(kwargs, :dispmode, "default"),
        "duty" => get(kwargs, :duty, ""),
        "dynadata" => get(kwargs, :dynadata, ""),
        "dynadll" => get(kwargs, :dynadll, "none"),
        "enabled" => get(kwargs, :enabled, true),
        "kv" => get(kwargs, :kv, 12.47),
        "kw" => get(kwargs, :kw, 0.0),
        "kva" => get(kwargs, :kva, 25.0),
        "kvar" => get(kwargs, :kvar, 0.0),
        "kwhrated" => get(kwargs, :kwhrated, 50.0),
        "kwhstored" => get(kwargs, :kwhstored, 50.0),
        "kwrated" => get(kwargs, :kwrated, 50.0),
        "model" => get(kwargs, :model, 1),
        "pf" => get(kwargs, :pf, 1.0),
        "phases" => get(kwargs, :phases, 3),
        "spectrum" => get(kwargs, :spectrum, "default"),
        "state" => get(kwargs, :state, "idling"),
        "timechargetrig" => get(kwargs, :timechargetrig, 2.0),
        "userdata" => get(kwargs, :userdata, ""),
        "usermodel" => get(kwargs, :usermodel, "none"),
        "vmaxpu" => get(kwargs, :vmaxpu, 1.1),
        "vminpu" => get(kwargs, :vimpu, 0.9),
        "yearly" => get(kwargs, :yearly, [1.0, 1.0]),
        "like" => get(kwargs, :like, "")
    )
end


"""
Creates a Dict{String,Any} containing all expected properties for a LoadShape
element. See OpenDSS documentation for valid fields and ways to specify
different properties.
"""
function _create_loadshape(name::String=""; kwargs...)
    if haskey(kwargs, :minterval)
        interval = kwargs[:minterval] / 60
    elseif haskey(kwargs, :sinterval)
        interval = kwargs[:sinterval] / 60 / 60
    else
        interval = get(kwargs, :interval, 1.0)
    end

    pmult = get(kwargs, :pmult, Vector{Float64}([]))
    qmult = get(kwargs, :qmult, pmult)

    npts = get(kwargs, :npts, length(pmult) == 0 && length(qmult) == 0 ? 0 : minimum(Int[length(a) for a in [pmult, qmult] if length(a) > 0]))

    pmult = pmult[1:npts]
    qmult = qmult[1:npts]

    hour = get(kwargs, :hour, collect(range(1.0, step=interval, length=npts)))[1:npts]

    Dict{String,Any}(
        "name" => name,
        "npts" => npts,
        "interval" => interval,
        "minterval" => interval * 60,
        "sinterval" => interval * 3600,
        "pmult" => pmult,
        "qmult" => qmult,
        "hour" => hour,
        "mean" => get(kwargs, :mean, mean(pmult)),
        "stddev" => get(kwargs, :stddev, std(pmult)),
        "csvfile" => get(kwargs, :csvfile, ""),
        "sngfile" => get(kwargs, :sngfile, ""),
        "dblfile" => get(kwargs, :dblfile, ""),
        "pqcsvfile" => get(kwargs, :pqcsvfile, ""),
        "action" => get(kwargs, :action, ""),
        "useactual" => get(kwargs, :useactual, true),
        "pmax" => get(kwargs, :pmax, 1.0),
        "qmax" => get(kwargs, :qmax, 1.0),
        "pbase" => get(kwargs, :pbase, 0.0),
        "like" => get(kwargs, :like, "")
    )
end


"""
Creates a Dict{String,Any} containing all expected properties for a XYCurve
object. See OpenDSS documentation for valid fields and ways to specify
different properties.
"""
function _create_xycurve(name::String=""; kwargs...)
    if haskey(kwargs, :points)
        xarray = Vector{Float64}([])
        yarray = Vector{Float64}([])

        i = 1
        for point in kwargs[:points]
            if i % 2 == 1
                push!(xarray, point)
            else
                push!(yarray, point)
            end
            i += 1
        end
    else
        xarray = get(kwargs, :xarray, Vector{Float64}([]))
        yarray = get(kwargs, :yarray, Vector{Float64}([]))
    end

    npts = min(length(xarray), length(yarray))

    points = Vector{Float64}([])
    for (x, y) in zip(xarray, yarray)
        push!(points, x)
        push!(points, y)
    end

    Dict{String,Any}(
        "name" => name,
        "npts" => npts,
        "points" => points,
        "yarray" => yarray,
        "xarray" => xarray,
        "csvfile" => get(kwargs, :csvfile, ""),
        "sngfile" => get(kwargs, :sngfile, ""),
        "dblfile" => get(kwargs, :dblfile, ""),
        "x" => get(kwargs, :x, isempty(xarray) ? 0.0 : xarray[1]),
        "y" => get(kwargs, :y, isempty(yarray) ? 0.0 : yarray[1]),
        "xshift" => get(kwargs, :xshift, 0),
        "yshift" => get(kwargs, :yshift, 0),
        "xscale" => get(kwargs, :xscale, 1.0),
        "yscale" => get(kwargs, :yscale, 1.0),
        "like" => get(kwargs, :like, ""),
    )
end


""
function _create_options(; kwargs...)
    Dict{String,Any}(
        "%growth" => get(kwargs, Symbol("%growth"), 2.5),
        "%mean" => get(kwargs, Symbol("%mean"), 65.0),
        "%normal" => get(kwargs, Symbol("%normal"), 100.0),
        "%stddev" => get(kwargs, Symbol("%stddev"), 9.0),
        "addtype" => get(kwargs, :addtype, "generator"),
        "algorithm" => get(kwargs, :algorithm, "newton"),
        "allocationfactors" => get(kwargs, :allocationfactors, ""),
        "allowduplicates" => get(kwargs, :allowduplicates, false),
        "autobuslist" => get(kwargs, :autobuslist, Vector{String}([])),
        "basefrequency" => get(kwargs, :basefrequency, 60.0),
        "bus" => get(kwargs, :bus, ""),
        "capkvar" => get(kwargs, :capkvar, 600.0),
        "casename" => get(kwargs, :casename, ""),
        "capmarkercode" => get(kwargs, :capmarkercode, 37),
        "capmarkersize" => get(kwargs, :capmarkersize, 3),
        "cfactors" => get(kwargs, :cfactors, 4.0),
        "circuit" => get(kwargs, :circuit, ""),
        "cktmodel" => get(kwargs, :cktmodel, "multiphase"),
        "class" => get(kwargs, :class, ""),
        "controlmode" => get(kwargs, :controlmode, "static"),
        "datapath" => get(kwargs, :datapath, ""),
        "defaultbasefrequency" => get(kwargs, :defaultbasefrequency, 60.0),
        "defaultbasefreq" => get(kwargs, :defaultbasefreq, 60.0), # Alias to defaultbasefrequency
        "defaultdaily" => get(kwargs, :defaultdaily, "default"),
        "defaultyearly" => get(kwargs, :defaultyearly, "default"),
        "demandinterval" => get(kwargs, :demandinterval, false),
        "diverbose" => get(kwargs, :diverbose, false),
        "dssvisualizationtool" => get(kwargs, :dssvisualizationtool, ""),
        "earthmodel" => get(kwargs, :earthmodel, "deri"),
        "editor" => get(kwargs, :editor, "notepad"),
        "element" => get(kwargs, :element, ""),
        "emergvmaxpu" => get(kwargs, :emergvmaxpu, 1.08),
        "emergvminpu" => get(kwargs, :emergvminpu, 0.90),
        "frequency" => get(kwargs, :frequency, 60.0),
        "genkw" => get(kwargs, :genkw, 1000.0),
        "genmult" => get(kwargs, :genmult, 1.0),
        "h" => get(kwargs, :h, ""),
        "harmonics" => get(kwargs, :harmonics, "all"),
        "hour" => get(kwargs, :hour, 1.0),
        "keeplist" => get(kwargs, :keeplist, Vector{String}([])),
        "ldcurve" => get(kwargs, :ldcurve, "nil"),
        "loadmodel" => get(kwargs, :loadmodel, "admittance"),
        "loadmult" => get(kwargs, :loadmult, 1.0),
        "log" => get(kwargs, :log, false),
        "lossregs" => get(kwargs, :lossregs, 13),
        "lossweight" => get(kwargs, :lossweight, 1.0),
        "markercode" => get(kwargs, :markercode, 0),
        "markswitches" => get(kwargs, :markswitches, false),
        "markcapacitors" => get(kwargs, :markcapacitors, false),
        "markpvsystems" => get(kwargs, :markpvsystems, false),
        "markregulators" => get(kwargs, :markregulators, false),
        "markstorage" => get(kwargs, :markstorage, false),
        "marktransformers" => get(kwargs, :marktransformers, false),
        "maxcontroliter" => get(kwargs, :maxcontroliter, 10),
        "maxiter" => get(kwargs, :maxiter, 15),
        "miniterations" => get(kwargs, :miniterations, 2),
        "mode" => get(kwargs, :mode, "Snap"),
        "name" => get(kwargs, :name, ""),
        "nodewidth" => get(kwargs, :nodewidth, 1),
        "normvmaxpu" => get(kwargs, :normvmaxpu, 1.05),
        "normvminpu" => get(kwargs, :normvminpu, 0.95),
        "numallociterations" => get(kwargs, :numallociterations, 2),
        "number" => get(kwargs, :number, 0),
        "object" => get(kwargs, :object, ""),
        "overloadreport" => get(kwargs, :overloadreport, false),
        "neglectloady" => get(kwargs, :neglectloady, false),
        "pricecurve" => get(kwargs, :pricecurve, ""),
        "pricesignal" => get(kwargs, :pricesignal, 25),
        "pvmarkercode" => get(kwargs, :pvmarkercode, 15),
        "pvmarkersize" => get(kwargs, :pvmarkersize, 1),
        "random" => get(kwargs, :random, "uniform"),
        "recorder" => get(kwargs, :recorder, false),
        "reduceoption" => get(kwargs, :reduceoption, "default"),
        "registryupdate" => get(kwargs, :registryupdate, true),
        "regmarkercode" => get(kwargs, :regmarkercode, 47),
        "regmarkersize" => get(kwargs, :regmarkersize, 1),
        "sampleenergymeters" => get(kwargs, :sampleenergymeters, false),
        "sec" => get(kwargs, :sec, 0.0),
        "showexport" => get(kwargs, :showexport, false),
        "stepsize" => get(kwargs, :stepsize, "1h"),
        "switchmarkercode" => get(kwargs, :switchmarkercode, 4),
        "terminal" => get(kwargs, :terminal, ""),
        "time" => get(kwargs, :time, Vector{Float64}([0.0, 0.0])),
        "tolerance" => get(kwargs, :tolerance, 0.0001),
        "totaltime" => get(kwargs, :totaltime, 0.0),
        "tracecontrol" => get(kwargs, :tracecontrol, false),
        "transmarkercode" => get(kwargs, :transmarkercode, 35),
        "transmarkersize" => get(kwargs, :transmarkersize, 1),
        "storemarkercode" => get(kwargs, :storemarkercode, 9),
        "storemarkersize" => get(kwargs, :storemarkersize, 1),
        "trapezoidal" => get(kwargs, :trapezoidal, false),
        "type" => get(kwargs, :type, ""),
        "ueregs" => get(kwargs, :ueregs, 11),
        "ueweight" => get(kwargs, :ueweight, 1.0),
        "voltagebases" => get(kwargs, :voltagebases, Vector{Float64}([])),
        "voltexceptionreport" => get(kwargs, :voltexceptionreport, false),
        "year" => get(kwargs, :year, 0),
        "zonelock" => get(kwargs, :zonelock, false),
    )
end



"Returns a Dict{String,Type} for the desired component `comp`, giving all of the expected data types"
const _dss_parameter_data_types = Dict{String,Dict{String,Type}}((comp, Dict{String,Type}((k, typeof(v)) for (k,v) in @eval $(Symbol("_create_$comp"))())) for comp in _dss_supported_components)
