"adds voltage balance indicators; should only be called after add_setpoint_bus_voltage!"
function add_setpoint_bus_voltage_balance_indicators!(pm::PMs.GenericPowerModel, sol)
    sol_dict = PMs.get(sol, "bus", Dict{String,Any}())

    num_conductors = length(PMs.conductor_ids(pm))
    @assert(PMs.ismulticonductor(pm) && num_conductors==3)

    if PMs.ismultinetwork(pm)
        bus_dict = pm.data["nw"]["$(pm.cnw)"]["bus"]
    else
        bus_dict = pm.data["bus"]
    end

    if length(bus_dict) > 0
        sol["bus"] = sol_dict
    end

    for (i,item) in bus_dict
        idx = Int(item["bus_i"])
        sol_item = sol_dict[i] = get(sol_dict, i, Dict{String,Any}())

        # assumes add_setpoint_bus_voltage! was called already
        vm = sol_item["vm"]
        va = sol_item["va"]
        v_abc = [vm[c]*exp(im*va[c]) for c in 1:3]
        a = exp(im*2*pi/3)
        v_pos = transpose(v_abc)*[1, a, a^2]./3
        v_neg = transpose(v_abc)*[1, a^2, a]./3
        v_zero = transpose(v_abc)*[1, 1, 1]./3

        sol_item["vm_seq_pos"] = abs(v_pos)
        sol_item["vm_seq_neg"] = abs(v_neg)
        sol_item["vm_seq_zero"] = abs(v_zero)
        sol_item["vuf"] = abs(v_neg)/abs(v_pos)
        sol_item["vm_ll"] = PMs.MultiConductorVector(abs.([1 -1 0; 0 1 -1; -1 0 1]*v_abc)./sqrt(3))
    end
end
