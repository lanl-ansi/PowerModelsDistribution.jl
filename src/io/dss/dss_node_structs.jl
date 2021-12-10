
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

    Dict{String,Any}(
        "name" => name,
        "phases" => get(kwargs, :phases, _get_implied_nphases(bus1; default=3)),
        "bus1" => bus1,
        "kv" => kv,
        "kw" => kw,
        "pf" => pf,
        "model" => get(kwargs, :model, 1),
        "yearly" => get(kwargs, :yearly, ""),
        "daily" => get(kwargs, :daily, ""),
        "duty" => get(kwargs, :duty, get(kwargs, :daily, "")),
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
        "zipv" => get(kwargs, :zipv, Float64[]),
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
    pf = get(kwargs, :pf, 0.80)

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
    elseif haskey(kwargs, :pf) && kwargs[:pf] != 0.80
        kvar = sign(pf) * kw * sqrt(1.0 / pf^2 - 1.0)
        kva = abs(kw) + kvar^2
    end

    kvarmax = get(kwargs, :maxkvar, kvar * 2.0)
    kvarmin = get(kwargs, :minkvar, -kvarmax)

    Dict{String,Any}(
        "name" => name,
        "phases" => get(kwargs, :phases, _get_implied_nphases(bus1; default=3)),
        "bus1" => bus1,
        "kv" => get(kwargs, :kv, 12.47),
        "kw" => kw,
        "pf" => pf,
        "model" => get(kwargs, :model, 1),
        "yearly" => get(kwargs, :yearly, get(kwargs, :daily, "")),
        "daily" => get(kwargs, :daily, ""),
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
PVSystem. See OpenDSS document
https://github.com/tshort/OpenDSS/blob/master/Doc/OpenDSS%20PVSystem%20Model.doc
for valid fields and ways to specify the different properties.
"""
function _create_pvsystem(name::String=""; kwargs...)
    bus1 = get(kwargs, :bus1, "")

    pmpp = get(kwargs, :pmpp, 500.0)
    if haskey(kwargs, Symbol("%pmpp"))
        pmpp *= kwargs[Symbol("%pmpp")] / 100.0
    end

    kv = get(kwargs, :kv, 12.47)
    pf = get(kwargs, :pf, 1.0)
    kvar = get(kwargs, :kvar, 0.0)
    kva = get(kwargs, :kva, 500.0)

    if haskey(kwargs, :pf) && !haskey(kwargs, :kva) && haskey(kwargs, :kvar)
        kva = kvar / sqrt(1 - pf^2) / sign(pf)
    elseif haskey(kwargs, :pf) && !haskey(kwargs, :kvar) && haskey(kwargs, :kva)
        kvar = kva * sqrt(1 - pf^2) * sign(pf)
    elseif !haskey(kwargs, :pf) && haskey(kwargs, :kva) && haskey(kwargs, :kvar)
        pf = sqrt(1 - (kvar / kva)^2) * sign(kvar)
    end

    Dict{String,Any}(
        "name" => name,
        "phases" => get(kwargs, :phases, _get_implied_nphases(bus1; default=3)),
        "bus1" => bus1,
        "kv" => kv,
        "irradiance" => get(kwargs, :irradiance, 1.0),
        "temperature" => get(kwargs, :temperature, 25.0),
        "pmpp" => get(kwargs, :pmpp, 500.0),
        "%pmpp" => get(kwargs, Symbol("%pmpp"), 100.0),
        "pf" => pf,
        "conn" => get(kwargs, :conn, WYE),
        "kvar" => kvar,
        "kva" => kva,
        "%cutin" => get(kwargs, Symbol("%cutin"), 20.0),
        "%cutout" => get(kwargs, Symbol("%cutout"), 20.0),
        "effcurve" => get(kwargs, :effcurve, ""),
        "p-tcurve" => get(kwargs, :ptcurve, ""),
        "%r" => get(kwargs, Symbol("%r"), 0),
        "%x" => get(kwargs, Symbol("%x"), 0.50),
        "model" => get(kwargs, :model, 1),
        "vminpu" => get(kwargs, :vminpu, 0.9),
        "vmaxpu" => get(kwargs, :vmaxpu, 1.1),
        "balanced" => get(kwargs, :balanced, false),
        "limitcurrent" => get(kwargs, :limitcurrent, false),
        "yearly" => get(kwargs, :yearly, get(kwargs, :daily, "")),
        "daily" => get(kwargs, :daily, ""),
        "duty" => get(kwargs, :duty, ""),
        "tyearly" => get(kwargs, :tyearly, ""),
        "tdaily" => get(kwargs, :tdaily, ""),
        "tduty" => get(kwargs, :tduty, ""),
        "class" => get(kwargs, :class, 1),
        "usermodel" => get(kwargs, :usermodel, ""),
        "userdata" => get(kwargs, :userdata, ""),
        "debugtrace" => get(kwargs, :debugtrace, false),
        "varfollowinverter" => get(kwargs, :varfollowinverter, false),
        "dutystart" => get(kwargs, :dutystart, 0),
        "wattpriority" => get(kwargs, :wattpriority, false),
        "pfpriority" => get(kwargs, :pfpriority, false),
        "%pminnovars" => get(kwargs, Symbol("%pminnovars"), -1.0),
        "%pminkvarmax" => get(kwargs, Symbol("%pminkvarmax"), -1.0),
        "kvarmax" => get(kwargs, :kvarmax, 500.0),
        "kvarmaxabs" => get(kwargs, :kvarmaxabs, 500.0),
        "spectrum" => get(kwargs, :spectrum, "defaultpvsystem"),
        "basefreq" => get(kwargs, :basefreq, 60.0),
        "enabled" => get(kwargs, :enabled, true),
        "like" => get(kwargs, :like, ""),
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
        "%charge" => get(kwargs, Symbol("%charge"), 100.0),
        "%discharge" => get(kwargs, Symbol("%discharge"), 100.0),
        "%effcharge" => get(kwargs, Symbol("%effcharge"), 90.0),
        "%effdischarge" => get(kwargs, Symbol("%effdischarge"), 90.0),
        "%idlingkvar" => get(kwargs, Symbol("%idlingkvar"), 0.0),
        "%idlingkw" => get(kwargs, Symbol("%idlingkw"), 1.0),
        "%r" => get(kwargs, Symbol("%r"), 0.0),
        "%reserve" => get(kwargs, Symbol("%reserve"), 20.0),
        "%stored" => get(kwargs, Symbol("%stored"), 100.0),
        "%x" => get(kwargs, Symbol("%x"), 50.0),
        "basefreq" => get(kwargs, :basefreq, 60.0),
        "bus1" => get(kwargs, :bus1, ""),
        "chargetrigger" => get(kwargs, :chargetrigger, 0.0),
        "class" => get(kwargs, :class, 0),
        "conn" => get(kwargs, :conn, WYE),
        "daily" => get(kwargs, :daily, ""),
        "debugtrace" => get(kwargs, :debugtrace, false),
        "dischargetrigger" => get(kwargs, :dischargetrigger, 0.0),
        "dispmode" => get(kwargs, :dispmode, "default"),
        "duty" => get(kwargs, :duty, ""),
        "dynadata" => get(kwargs, :dynadata, ""),
        "dynadll" => get(kwargs, :dynadll, "none"),
        "enabled" => get(kwargs, :enabled, true),
        "kv" => get(kwargs, :kv, 12.47),
        "kw" => get(kwargs, :kw, 0.0),
        "kva" => get(kwargs, :kva, get(kwargs, :kwrated, 50.0)),
        "kvar" => get(kwargs, :kvar, sign(get(kwargs, :pf, 1.0)) * get(kwargs, :kw, 0.0) * sqrt(1.0 / get(kwargs, :pf, 1.0)^2 - 1.0)),
        "kwhrated" => get(kwargs, :kwhrated, 50.0),
        "kwhstored" => get(kwargs, :kwhstored, get(kwargs, :kwhrated, 50.0) * get(kwargs, :stored, 100.0) / 100.0),
        "kwrated" => get(kwargs, :kwrated, 50.0),
        "model" => get(kwargs, :model, 1),
        "pf" => get(kwargs, :pf, get(kwargs, :kw, 0.0) == 0 ? 1.0 : get(kwargs, :kvar, kwargs[:kw]) / kwargs[:kw]),
        "phases" => get(kwargs, :phases, _get_implied_nphases(get(kwargs, :bus1, ""); default=3)),
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
