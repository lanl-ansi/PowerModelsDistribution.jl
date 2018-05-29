"variable: `w[i] >= 0` for `i` in `bus`es"
function variable_tp_voltage_magnitude_sqr(pm::GenericPowerModel; nw::Int=pm.cnw, ph::Int=pm.cph, bounded=true)
    bus_ph = [(i, h) for i in ids(pm, nw, :bus) for h in PMs.phase_ids(pm)]

    if bounded
        var(pm, nw)[:w] = @variable(pm.model,
            [i in bus_ph], basename="$(nw)_$(ph)_w",
            lowerbound = ref(pm, nw, :bus, i[1], "vmin", i[2])^2,
            upperbound = ref(pm, nw, :bus, i[1], "vmax", i[2])^2,
            start = PMs.getval(ref(pm, nw, :bus, i[1]), "w_start", i[2], 1.001)
        )
    else
        var(pm, nw)[:w] = @variable(pm.model,
            [i in bus_ph], basename="$(nw)_$(ph)_w",
            lowerbound = 0,
            start = PMs.getval(ref(pm, nw, :bus, i[1]), "w_start", i[2], 1.001)
        )
    end
end


""
function variable_tp_voltage_product(pm::GenericPowerModel; nw::Int=pm.cnw, ph::Int=pm.cph, bounded=true)
    if bounded
        wr_min, wr_max, wi_min, wi_max = PMs.calc_voltage_product_bounds(ref(pm, nw, :buspairs), ph)
        bp_phf_pht = [(i..., h, g) for i in keys(ref(pm, nw, :buspairs)) for h in PMs.phase_ids(pm) for g in PMs.phase_ids(pm)]
        bus_ph = [(i, i, h, g) for i in ids(pm, nw, :bus) for h in PMs.phase_ids(pm) for g in PMs.phase_ids(pm) if h != g]

        append!(bus_ph, bp_phf_pht)

        var(pm, nw)[:wr] = @variable(pm.model,
            [b in bus_ph], basename="$(nw)_wr",
            lowerbound = b[1] != b[2] ? wr_min[b[1:2]] : ref(pm, nw, :bus, b[1], "vmin", b[3])^2,  # TODO: correct bounds for different phases
            upperbound = b[1] != b[2] ? wr_max[b[1:2]] : ref(pm, nw, :bus, b[1], "vmax", b[3])^2,  # TODO: correct bounds for different phases
            start = PMs.getval(b[1] != b[2] ? ref(pm, nw, :buspairs, b[1:2]) : ref(pm, nw, :bus, b[1]), "wr_start", b[3], 1.0)
        )
        var(pm, nw)[:wi] = @variable(pm.model,
            [b in bus_ph], basename="$(nw)_wi",
            lowerbound = b[1] != b[2] ? wi_min[b[1:2]] : ref(pm, nw, :bus, b[1], "vmin", b[3])^2,  # TODO: correct bounds for different phases
            upperbound = b[1] != b[2] ? wi_max[b[1:2]] : ref(pm, nw, :bus, b[1], "vmax", b[3])^2,  # TODO: correct bounds for different phases
            start = PMs.getval(b[1] != b[2] ? ref(pm, nw, :buspairs, b[1:2]) : ref(pm, nw, :bus, b[1]), "wi_start", b[3])
        )
    else
        bp_phf_pht = [(i..., h, g) for i in keys(ref(pm, nw, :buspairs)) for h in PMs.phase_ids(pm) for g in PMs.phase_ids(pm)]
        bus_ph = [(i, i, h, g) for i in ids(pm, nw, :bus) for h in PMs.phase_ids(pm) for g in PMs.phase_ids(pm) if g != h]

        append!(bus_ph, bp_phf_pht)

        var(pm, nw)[:wr] = @variable(pm.model,
            [b in bus_ph], basename="$(nw)_wr",
            start = PMs.getval(b[1] != b[2] ? ref(pm, nw, :buspairs, b[1:2]) : ref(pm, nw, :bus, b[1]), "wr_start", b[3], 1.0)
        )
        var(pm, nw)[:wi] = @variable(pm.model,
            [b in bus_ph], basename="$(nw)_wi",
            start = PMs.getval(b[1] != b[2] ? ref(pm, nw, :buspairs, b[1:2]) : ref(pm, nw, :bus, b[1]), "wi_start", b[3])
        )
    end
end