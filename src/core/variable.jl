"variable: `w[i] >= 0` for `i` in `bus`es"
function variable_tp_voltage_magnitude_sqr(pm::GenericPowerModel; nw::Int=pm.cnw, ph::Int=pm.cph, bounded=true)
    bus_ph = [(i, h) for i in ids(pm, nw, :bus) for h in PMs.phase_ids(pm)]

    if bounded
        W = var(pm, nw)[:w] = @variable(pm.model,
            [i in bus_ph], basename="$(nw)_w",
            lowerbound = ref(pm, nw, :bus, i[1], "vmin", i[2])^2,
            upperbound = ref(pm, nw, :bus, i[1], "vmax", i[2])^2,
            start = PMs.getval(ref(pm, nw, :bus, i[1]), "w_start", i[2], 1.001)
        )
    else
        W = var(pm, nw)[:w] = @variable(pm.model,
            [i in bus_ph], basename="$(nw)_w",
            lowerbound = 0,
            start = PMs.getval(ref(pm, nw, :bus, i[1]), "w_start", i[2], 1.001)
        )
    end

    var(pm, nw, ph)[:w] = Dict{Int,Any}()
    for i in ids(pm, nw, :bus)
        var(pm, nw, ph, :w)[i] = W[(i, ph)]
    end
end


""
function variable_tp_voltage_product(pm::GenericPowerModel; nw::Int=pm.cnw, ph::Int=pm.cph, bounded=true)
    bp_phf_pht = [(i, j, h, g) for (i,j) in keys(ref(pm, nw, :buspairs)) for h in PMs.phase_ids(pm) for g in PMs.phase_ids(pm)]
    bus_ph = [(i, i, h, g) for i in ids(pm, nw, :bus) for h in PMs.phase_ids(pm) for g in PMs.phase_ids(pm) if h != g]
    append!(bus_ph, bp_phf_pht)

    WR = var(pm, nw)[:wr] = @variable(pm.model,
        [b in bus_ph], basename="$(nw)_wr",
        start = PMs.getval(b[1] != b[2] ? ref(pm, nw, :buspairs, b[1:2]) : ref(pm, nw, :bus, b[1]), "wr_start", b[3], 1.0)
    )

    WI = var(pm, nw)[:wi] = @variable(pm.model,
        [b in bus_ph], basename="$(nw)_wi",
        start = PMs.getval(b[1] != b[2] ? ref(pm, nw, :buspairs, b[1:2]) : ref(pm, nw, :bus, b[1]), "wi_start", b[3])
    )

    if bounded
        # Diagonal bounds
        wr_min, wr_max, wi_min, wi_max = PMs.calc_voltage_product_bounds(ref(pm, nw, :buspairs), ph)
        for (i, j) in ids(pm, nw, :buspairs)
            PMs.setupperbound(WR[(i, j, ph, ph)], wr_max[(i,j)])
            PMs.setupperbound(WI[(i, j, ph, ph)], wi_max[(i,j)])

            PMs.setlowerbound(WR[(i, j, ph, ph)], wr_min[(i,j)])
            PMs.setlowerbound(WI[(i, j, ph, ph)], wi_min[(i,j)])
        end

        # Off-diagonal bounds
        for h in PMs.phase_ids(pm)
            if h != ph
                wr_min, wr_max, wi_min, wi_max = calc_tp_voltage_product_bounds(pm, bus_ph)
                for k in bus_ph
                    PMs.setupperbound(WR[k], wr_max[k])
                    PMs.setupperbound(WI[k], wi_max[k])

                    PMs.setlowerbound(WR[k], wr_min[k])
                    PMs.setlowerbound(WI[k], wi_min[k])
                end
            end
        end
    end

    var(pm, nw, ph)[:wr] = Dict{Tuple{Int,Int},Any}()
    var(pm, nw, ph)[:wi] = Dict{Tuple{Int,Int},Any}()
    for (i, j) in ids(pm, nw, :buspairs)
        var(pm, nw, ph, :wr)[(i,j)] = WR[(i, j, ph, ph)]
        var(pm, nw, ph, :wi)[(i,j)] = WI[(i, j, ph, ph)]
    end
end
