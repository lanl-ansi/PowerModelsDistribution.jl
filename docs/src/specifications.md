# Problem Specifications


## Optimal Power Flow (OPF) with On-Load Tap Changers (OLTC)
This problem is identical to `tp_opf`, except that all transformers are now modelled as on-load tap changers (OLTCs). Each phase has an individual tap ratio, which can be either variable or fixed, as specified in the data model.
### Objective
```julia
objective_min_fuel_cost(pm)
```

### Variables
```julia
variable_tp_voltage(pm)
variable_tp_branch_flow(pm)

for c in PMs.conductor_ids(pm)
    PMs.variable_generation(pm, cnd=c)
    PMs.variable_dcline_flow(pm, cnd=c)
end
variable_tp_trans_flow(pm)
variable_tp_oltc_tap(pm)
```

### Constraints
```julia
constraint_tp_voltage(pm)

for i in PMs.ids(pm, :ref_buses)
    constraint_tp_theta_ref(pm, i)
end

for i in PMs.ids(pm, :bus), c in PMs.conductor_ids(pm)
    constraint_kcl_shunt_trans(pm, i, cnd=c)
end

for i in PMs.ids(pm, :branch)
    for c in PMs.conductor_ids(pm)
        constraint_ohms_tp_yt_from(pm, i, cnd=c)
        constraint_ohms_tp_yt_to(pm, i, cnd=c)

        PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

        PMs.constraint_thermal_limit_from(pm, i, cnd=c)
        PMs.constraint_thermal_limit_to(pm, i, cnd=c)
    end
end

for i in PMs.ids(pm, :dcline), c in PMs.conductor_ids(pm)
    PMs.constraint_dcline(pm, i, cnd=c)
end

for i in PMs.ids(pm, :trans)
    constraint_tp_oltc(pm, i)
end
```

## Optimal Power Flow (OPF) with Load Models (LM)
Unlike `tp_opf`, which models all loads as constant power loads, this problem specification additionally supports loads proportional to the voltage magnitude (a.k.a. constant current) and the square of the voltage magnitude (a.k.a. constant impedance). Each load now has associated active and reactive power variables. In `tp_opf`, loads are directly added as parameters in KCL.

### Objective
```julia
objective_min_fuel_cost(pm)
```

### Variables
```julia
variable_tp_voltage(pm)
variable_tp_branch_flow(pm)

for c in PMs.conductor_ids(pm)
    PMs.variable_generation(pm, cnd=c)
    PMs.variable_dcline_flow(pm, cnd=c)
end
variable_tp_trans_flow(pm)
variable_tp_oltc_tap(pm)
```

### Constraints
```julia
constraint_tp_voltage(pm)

for i in PMs.ids(pm, :ref_buses)
    constraint_tp_theta_ref(pm, i)
end

for i in PMs.ids(pm, :bus), c in PMs.conductor_ids(pm)
    constraint_kcl_shunt_trans_load(pm, i, cnd=c)
end

for id in PMs.ids(pm, :load)
    model = PMs.ref(pm, pm.cnw, :load, id, "model")
    if model=="constant_power"
        constraint_tp_load_power_setpoint(pm, id)
    elseif model=="proportional_vm"
        constraint_tp_load_power_prop_vm(pm, id)
    elseif model=="proportional_vmsqr"
        constraint_tp_load_power_prop_vmsqr(pm, id)
    else
        Memento.@error(LOGGER, "Unknown model $model for load $id.")
    end
end

for i in PMs.ids(pm, :branch)
    for c in PMs.conductor_ids(pm)
        constraint_ohms_tp_yt_from(pm, i, cnd=c)
        constraint_ohms_tp_yt_to(pm, i, cnd=c)

        PMs.constraint_voltage_angle_difference(pm, i, cnd=c)

        PMs.constraint_thermal_limit_from(pm, i, cnd=c)
        PMs.constraint_thermal_limit_to(pm, i, cnd=c)
    end
end

for i in PMs.ids(pm, :dcline), c in PMs.conductor_ids(pm)
    PMs.constraint_dcline(pm, i, cnd=c)
end

for i in PMs.ids(pm, :trans)
    constraint_tp_trans(pm, i)
end
```

## Power Flow (PF) with Load Models (LM)
Unlike `tp_pf`, which models all loads as constant power loads, this problem specification additionally supports loads proportional to the voltage magnitude (a.k.a. constant current) and the square of the voltage magnitude (a.k.a. constant impedance). Each load now has associated active and reactive power variables. In `tp_pf`, loads are directly added as parameters in KCL.

### Variables
```julia
variable_tp_voltage(pm, bounded=false)
variable_tp_branch_flow(pm, bounded=false)

for c in PMs.conductor_ids(pm)
    PMs.variable_generation(pm, bounded=false, cnd=c)
    variable_load_flow(pm, cnd=c)
    PMs.variable_dcline_flow(pm, bounded=false, cnd=c)
end

variable_tp_trans_flow(pm, bounded=false)
```

### Constraints
```julia
constraint_tp_voltage(pm, bounded=false)

for (i,bus) in PMs.ref(pm, :ref_buses)
    constraint_tp_theta_ref(pm, i)

    for c in PMs.conductor_ids(pm)
        @assert bus["bus_type"] == 3
        PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c)
    end
end

for (i,bus) in PMs.ref(pm, :bus), c in PMs.conductor_ids(pm)
    constraint_kcl_shunt_trans_load(pm, i, cnd=c)

    # PV Bus Constraints
    if length(PMs.ref(pm, :bus_gens, i)) > 0 && !(i in PMs.ids(pm,:ref_buses))
        # this assumes inactive generators are filtered out of bus_gens
        @assert bus["bus_type"] == 2

        PMs.constraint_voltage_magnitude_setpoint(pm, i, cnd=c)
        for j in PMs.ref(pm, :bus_gens, i)
            PMs.constraint_active_gen_setpoint(pm, j, cnd=c)
        end
    end
end

for id in PMs.ids(pm, :load)
    model = PMs.ref(pm, pm.cnw, :load, id, "model")
    if model=="constant_power"
        constraint_tp_load_power_setpoint(pm, id)
    elseif model=="proportional_vm"
        constraint_tp_load_power_prop_vm(pm, id)
    elseif model=="proportional_vmsqr"
        constraint_tp_load_power_prop_vmsqr(pm, id)
    else
        Memento.@error(LOGGER, "Unknown model $model for load $id.")
    end
end

for i in PMs.ids(pm, :branch), c in PMs.conductor_ids(pm)
    constraint_ohms_tp_yt_from(pm, i, cnd=c)
    constraint_ohms_tp_yt_to(pm, i, cnd=c)
    # PMs.constraint_ohms_yt_from(pm, i, cnd=c)
    # PMs.constraint_ohms_yt_to(pm, i, cnd=c)
end

for (i,dcline) in PMs.ref(pm, :dcline), c in PMs.conductor_ids(pm)

    PMs.constraint_active_dcline_setpoint(pm, i, cnd=c)

    f_bus = PMs.ref(pm, :bus)[dcline["f_bus"]]
    if f_bus["bus_type"] == 1
        PMs.constraint_voltage_magnitude_setpoint(pm, f_bus["index"], cnd=c)
    end

    t_bus = PMs.ref(pm, :bus)[dcline["t_bus"]]
    if t_bus["bus_type"] == 1
        PMs.constraint_voltage_magnitude_setpoint(pm, t_bus["index"], cnd=c)
    end
end

for i in PMs.ids(pm, :trans)
    constraint_tp_trans(pm, i)
end
```
