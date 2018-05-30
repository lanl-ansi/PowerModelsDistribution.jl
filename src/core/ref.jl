""
function calc_tp_voltage_product_bounds(pm::GenericPowerModel, buspairs; nw::Int=pm.cnw)
    wr_min = Dict([(bp, -Inf) for bp in buspairs])
    wr_max = Dict([(bp,  Inf) for bp in buspairs])
    wi_min = Dict([(bp, -Inf) for bp in buspairs])
    wi_max = Dict([(bp,  Inf) for bp in buspairs])

    for (i, j, h, g) in buspairs
        if i == j
            bus = ref(pm, nw, :bus)[i]
            angmin = -1.0472
            angmax =  1.0472
            vm_fr_max = getmpv(bus["vmax"], h)
            vm_to_max = getmpv(bus["vmax"], g)
            vm_fr_min = getmpv(bus["vmin"], h)
            vm_to_min = getmpv(bus["vmin"], g)
        else
            buspair = ref(pm, nw, :buspairs)[(i, j)]
            angmin = min(getmpv(buspair["angmin"], h), getmpv(buspair["angmin"], g))
            angmax = max(getmpv(buspair["angmax"], h), getmpv(buspair["angmax"], g))
            vm_fr_max = getmpv(buspair["vm_fr_max"], h)
            vm_to_max = getmpv(buspair["vm_to_max"], g)
            vm_fr_min = getmpv(buspair["vm_fr_min"], h)
            vm_to_min = getmpv(buspair["vm_to_min"], g)
        end

        if angmin >= 0
            wr_max[(i, j, h, g)] = vm_fr_max * vm_to_max * cos(angmin)
            wr_min[(i, j, h, g)] = vm_fr_min * vm_to_min * cos(angmax)
            wi_max[(i, j, h, g)] = vm_fr_max * vm_to_max * sin(angmax)
            wi_min[(i, j, h, g)] = vm_fr_min * vm_to_min * sin(angmin)
        end
        if angmax <= 0
            wr_max[(i, j, h, g)] = vm_fr_max * vm_to_max * cos(angmax)
            wr_min[(i, j, h, g)] = vm_fr_min * vm_to_min * cos(angmin)
            wi_max[(i, j, h, g)] = vm_fr_min * vm_to_min * sin(angmax)
            wi_min[(i, j, h, g)] = vm_fr_max * vm_to_max * sin(angmin)
        end
        if angmin < 0 && angmax > 0
            wr_max[(i, j, h, g)] = vm_fr_max * vm_to_max * 1.0
            wr_min[(i, j, h, g)] = vm_fr_min * vm_to_min * min(cos(angmin), cos(angmax))
            wi_max[(i, j, h, g)] = vm_fr_max * vm_to_max * sin(angmax)
            wi_min[(i, j, h, g)] = vm_fr_max * vm_to_max * sin(angmin)
        end

    end

    return wr_min, wr_max, wi_min, wi_max
end