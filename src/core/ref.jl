""
function calc_tp_voltage_product_bounds(pm::GenericPowerModel, buspairs; nw::Int=pm.cnw)
    wr_min = Dict([(bp, -Inf) for bp in buspairs])
    wr_max = Dict([(bp,  Inf) for bp in buspairs])
    wi_min = Dict([(bp, -Inf) for bp in buspairs])
    wi_max = Dict([(bp,  Inf) for bp in buspairs])

    for (i, j, h, g) in buspairs
        if i == j
            bus = ref(pm, nw, :bus)[i]
            vm_fr_max = bus["vmax"][h]
            vm_to_max = bus["vmax"][g]
            vm_fr_min = bus["vmin"][h]
            vm_to_min = bus["vmin"][g]
        else
            buspair = ref(pm, nw, :buspairs)[(i, j)]
            vm_fr_max = buspair["vm_fr_max"][h]
            vm_to_max = buspair["vm_to_max"][g]
            vm_fr_min = buspair["vm_fr_min"][h]
            vm_to_min = buspair["vm_to_min"][g]
        end

        wr_max[(i, j, h, g)] =  vm_fr_max * vm_to_max
        wr_min[(i, j, h, g)] = -vm_fr_max * vm_to_max
        wi_max[(i, j, h, g)] =  vm_fr_max * vm_to_max
        wi_min[(i, j, h, g)] = -vm_fr_max * vm_to_max
    end

    return wr_min, wr_max, wi_min, wi_max
end