function all_gens_on(result)
    # tolerance of 1e-5 is needed for SCS tests to pass
    return minimum([gen["gen_status"] for (i,gen) in result["solution"]["gen"]]) >= 1.0 - 1e-5
end

function active_power_served(result)
    return sum([load["pd"] for (i,load) in result["solution"]["load"]])
end

function all_voltages_on(result)
    return return minimum([bus["status"] for (i,bus) in result["solution"]["bus"]]) >= 1.0 - 1e-3 #(note, non-SCS solvers are more accurate)
end


# helper functions to access solutions by their OpenDSS names
bus_name2id(pmd_data, name) = [bus["index"] for (_,bus) in pmd_data["bus"] if haskey(bus, "name") && bus["name"]==name][1]
va(sol, pmd_data, name) = PMD._wrap_to_pi(sol["solution"]["bus"][string(bus_name2id(pmd_data, name))]["va"][:])*180/pi
vm(sol, pmd_data, name) = sol["solution"]["bus"][string(bus_name2id(pmd_data, name))]["vm"]
tap(i, pm) = JuMP.value.(PMD.var(pm, InfrastructureModels.nw_id_default, :tap)[i])
vi(sol, pmd_data, name) = sol["solution"]["bus"][string(bus_name2id(pmd_data, name))]["vi"]
vr(sol, pmd_data, name) = sol["solution"]["bus"][string(bus_name2id(pmd_data, name))]["vr"]
calc_vm_acr(sol, pmd_data, name) = sqrt.(vi(sol, pmd_data, name).^2 .+ vr(sol, pmd_data, name).^2)
calc_va_acr(sol, pmd_data, name) = rad2deg.(PMD._wrap_to_pi(atan.(vi(sol, pmd_data, name), vr(sol, pmd_data, name))))

# Helper functions adjusted for eng model
vi(sol, name) = sol["solution"]["bus"][name]["vi"]
vr(sol, name) = sol["solution"]["bus"][name]["vr"]
calc_vm_acr(sol, name) = sqrt.(vi(sol, name).^2 .+ vr(sol, name).^2)
calc_va_acr(sol, name) = rad2deg.(PMD._wrap_to_pi(atan.(vi(sol, name), vr(sol, name))))
va(sol, name) = PMD._wrap_to_pi(sol["solution"]["bus"][name]["va"][:])*180/pi
vm(sol, name) = sol["solution"]["bus"][name]["vm"]
pd(sol, name) = sol["solution"]["load"][name]["pd_bus"]
qd(sol, name) = sol["solution"]["load"][name]["qd_bus"]


# Helper functions for load models tests
load_name2id(pmd_data, name) = [load["index"] for (_,load) in pmd_data["load"] if haskey(load, "name") && load["name"]==name][1]
pdvar(pm, pmd_data, name) = [PMD.var(pm, nw_id_default, c, :pd, load_name2id(pmd_data, name)) for c in 1:3]
pd(sol, pmd_data, name) = sol["solution"]["load"][string(load_name2id(pmd_data, name))]["pd_bus"]
qdvar(pm, pmd_data, name) = [PMD.var(pm, nw_id_default, c, :qd, load_name2id(pmd_data, name)) for c in 1:3]
qd(sol, pmd_data, name) = sol["solution"]["load"][string(load_name2id(pmd_data, name))]["qd_bus"]
sd(pm, pmd_data, name) = pd(sol, pmd_data, name)+im*qd(sol, pmd_data, name)

calc_vm_w(result, id) = sqrt.(      result["solution"]["bus"][id]["w"])
calc_vm_W(result, id) = sqrt.(LinearAlgebra.diag( result["solution"]["bus"][id]["Wr"]))


"The solar parsing is a bit off, this method corrects for that."
function pv1_correction!(data_eng)
    pv1 = data_eng["solar"]["pv1"]
    pv1["pg_lb"] = pv1["pg_ub"] = pv1["pg"]
    pv1["qg_lb"] = pv1["qg_ub"] = pv1["qg"]
end


"""
The ACR formulation returns non-physical solutions when the neutral voltage is not bounded below.
This method uses the ODD solution to find valid lower bounds, so the formulation can be validated.
"""
function add_neutral_lb_from_soldss(data_eng, sol_dss)
    de = deepcopy(data_eng)
    for (id, sol_bus) in sol_dss["bus"]
        eng_bus = de["bus"][id]
        ts = eng_bus["terminals"]
        for t in [t for t in ts if t>3 && t in keys(sol_bus["vm"])]
            idx = findfirst(ts.==t)
            eng_bus["vm_lb"] = haskey(eng_bus, "vm_lb") ? eng_bus["vm_lb"] : fill(0.0, length(ts))
            eng_bus["vm_lb"][idx] = max(eng_bus["vm_lb"][idx], sol_bus["vm"][t]*1E-3*0.9)
        end
    end
    return de
end


"Compares a PMD and OpenDSS solution, and returns the largest difference in voltage profile in per unit."
function compare_sol_dss_pmd(sol_dss::Dict{String,<:Any}, sol_pmd::Dict{String,<:Any}, data_eng::Dict{String,<:Any}, data_math::Dict{String,<:Any}; compare_math::Bool=false, verbose::Bool=true, floating_buses::Vector=[], skip_buses::Vector=[], v_err_print_tol::Real=1E-6)
    max_v_err_pu = 0.0

    # voltage base for ENGINEERING buses in [V]
    vbase = Dict(id=>data_math["bus"]["$ind"]["vbase"]*data_math["settings"]["voltage_scale_factor"] for (id,ind) in data_math["bus_lookup"])

    buses_intersected = intersect(keys(sol_dss["bus"]), keys(sol_pmd["bus"]))
    for id in setdiff(buses_intersected, skip_buses)
        pmd_bus = sol_pmd["bus"][id]
        dss_bus = sol_dss["bus"][id]

        terminals = data_eng["bus"][id]["terminals"]
        if compare_math
            ts = filter(x->haskey(dss_bus["vm"], "$x"), terminals)
            v_dss = [dss_bus["vm"]["$t"]*exp(im*dss_bus["va"]["$t"]) for t in ts]
            # convert to V instead of usual kV
            v_pmd = [pmd_bus["vm"][idx]*exp(im*deg2rad(pmd_bus["va"][idx]))*data_eng["settings"]["voltage_scale_factor"] for (idx,t) in enumerate(ts)]
        else
            ts = filter(x->haskey(dss_bus["vm"], x), terminals)
            v_dss = [dss_bus["vm"]["$t"]*exp(im*dss_bus["va"]["$t"]) for t in ts]
            # convert to V instead of usual kV
            v_pmd = [(pmd_bus["vr"][idx]+im*pmd_bus["vi"][idx])*data_eng["settings"]["voltage_scale_factor"] for (idx,t) in enumerate(ts)]
        end

        # convert to pu
        v_dss_pu = v_dss/vbase[id]
        v_pmd_pu = v_pmd/vbase[id]

        # convert to diffs if floating
        N = length(v_dss)
        if id in floating_buses && N>1
            v_dss_pu = [v_dss_pu[i]-v_dss_pu[j] for i in 1:N for j in i:N if i!=j]/vbase[id]
            v_pmd_pu = [v_pmd_pu[i]-v_pmd_pu[j] for i in 1:N for j in i:N if i!=j]/vbase[id]
            labels = ["$(ts[i])-$(ts[j])" for i in 1:N for j in i:N if i!=j]
        else
            labels = string.(ts)
        end

        for i in eachindex(v_pmd_pu)
            v_err_pu = abs.(v_dss_pu[i]-v_pmd_pu[i]); max_v_err_pu = max(max_v_err_pu, v_err_pu)

            if v_err_pu>v_err_print_tol && verbose
                println("terminal $id.$(labels[i])")
                println("\t |U| dss: $(abs(v_dss_pu[i]))")
                println("\t     pmd: $(abs(v_pmd_pu[i]))")
                println("\t  ∠U dss:  $(angle(v_dss_pu[i]))")
                println("\t     pmd:  $(angle(v_pmd_pu[i]))")
            end
        end
    end

    return max_v_err_pu
end
