"adds voltage balance indicators; should only be called after add_setpoint_bus_voltage!"
function add_setpoint_bus_voltage_balance_indicators!(pm::_PMs.AbstractPowerModel, sol)
    sol_dict = _PMs.get(sol, "bus", Dict{String,Any}())

    num_conductors = length(_PMs.conductor_ids(pm))
    @assert(_PMs.ismulticonductor(pm) && num_conductors==3)

    if _PMs.ismultinetwork(pm)
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
        sol_item["vm_ll"] = _PMs.MultiConductorVector(abs.([1 -1 0; 0 1 -1; -1 0 1]*v_abc)./sqrt(3))
    end
end


""
function solution_pbs!(pm::_PMs.AbstractPowerModel, sol::Dict{String,Any})
    _PMs.add_setpoint_bus_voltage!(sol, pm)
    _PMs.add_setpoint_generator_power!(sol, pm)
    _PMs.add_setpoint_branch_flow!(sol, pm)
    add_setpoint_bus_slack!(sol, pm)
end


""
function add_setpoint_bus_slack!(sol, pm::_PMs.AbstractPowerModel)
    _PMs.add_setpoint!(sol, pm, "bus", "p_slack", :p_slack, status_name="bus_type", inactive_status_value=4)
    _PMs.add_setpoint!(sol, pm, "bus", "q_slack", :q_slack, status_name="bus_type", inactive_status_value=4)
end


""
function solution_bctr!(pm::_PMs.AbstractPowerModel, sol::Dict{String,<:Any})
    _PMs.add_setpoint_bus_voltage!(sol, pm)
    add_setpoint_bus_voltage_balance_indicators!(pm, sol)
    _PMs.add_setpoint_generator_power!(sol, pm)
end


""
function solution_tp!(pm::_PMs.AbstractPowerModel, sol::Dict{String,Any})
    add_setpoint_bus_voltage!(sol, pm)
    _PMs.add_setpoint_generator_power!(sol, pm)
    add_setpoint_branch_flow!(sol, pm)
    add_setpoint_branch_current!(sol, pm)
    _PMs.add_setpoint_dcline_flow!(sol, pm)

    _PMs.add_dual_kcl!(sol, pm)
    _PMs.add_dual_sm!(sol, pm) # Adds the duals of the transmission lines' thermal limits.

    if haskey(pm.setting, "output") && haskey(pm.setting["output"], "original_variables") && pm.setting["output"]["original_variables"] == true
        add_rank!(sol, pm)
        add_is_ac_feasible!(sol, pm)
        add_original_variables!(sol, pm)
    end
end


""
function add_setpoint_bus_voltage!(sol, pm::_PMs.AbstractPowerModel)
    _PMs.add_setpoint!(sol, pm, "bus", "vm", :w; scale = (x,item,i) -> sqrt(x), status_name="bus_type", inactive_status_value=4)
    _PMs.add_setpoint!(sol, pm, "bus", "w",  :w, status_name="bus_type", inactive_status_value=4)
    _PMs.add_setpoint!(sol, pm, "bus", "wr", :wr, status_name="bus_type", inactive_status_value=4)
    _PMs.add_setpoint!(sol, pm, "bus", "wi", :wi, status_name="bus_type", inactive_status_value=4)
end


""
function add_setpoint_branch_flow!(sol, pm::_PMs.AbstractPowerModel)
    if haskey(pm.setting, "output") && haskey(pm.setting["output"], "branch_flows") && pm.setting["output"]["branch_flows"] == true
        _PMs.add_setpoint!(sol, pm, "branch", "pf", :p; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])], status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "qf", :q; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])], status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "pf_ut", :p_ut; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])], status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "qf_ut", :q_ut; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])], status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "pf_lt", :p_lt; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])], status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "qf_lt", :q_lt; extract_var = (var,idx,item) -> var[(idx, item["f_bus"], item["t_bus"])], status_name="br_status")

        _PMs.add_setpoint!(sol, pm, "branch", "pt", :p; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])], status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "qt", :q; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])], status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "pt_ut", :p_ut; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])], status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "qt_ut", :q_ut; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])], status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "pt_lt", :p_lt; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])], status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "qt_lt", :q_lt; extract_var = (var,idx,item) -> var[(idx, item["t_bus"], item["f_bus"])], status_name="br_status")

    end
end


""
function add_setpoint_branch_current!(sol, pm::_PMs.AbstractPowerModel)
    if haskey(pm.setting, "output") && haskey(pm.setting["output"], "branch_flows") && pm.setting["output"]["branch_flows"] == true
        _PMs.add_setpoint!(sol, pm, "branch", "ccm", :ccm; scale = (x,item,i) -> sqrt(x), status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "cc", :ccm, status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "ccr", :ccmr, status_name="br_status")
        _PMs.add_setpoint!(sol, pm, "branch", "cci", :ccmi, status_name="br_status")
    end
end


""
function add_rank!(sol, pm::_PMs.AbstractPowerModel; tol = 1e-6)
    add_rank_voltage_variable!(sol, pm; tol=tol)
    add_rank_current_variable!(sol, pm; tol=tol)
    add_rank_branch_flow!(sol, pm; tol=tol)
end


""
function add_rank_voltage_variable!(sol, pm::_PMs.AbstractPowerModel; tol = 1e-6)
end


""
function add_rank_current_variable!(sol, pm::_PMs.AbstractPowerModel; tol = 1e-6)
end


""
function add_rank_branch_flow!(sol, pm::_PMs.AbstractPowerModel; tol = 1e-6)
end


""
function add_is_ac_feasible!(sol, pm::_PMs.AbstractPowerModel)
end


""
function add_is_ac_feasible!(sol, pm::AbstractUBFModels)
    for (nw, network) in pm._PMs.ref[:nw]
        branch_feasibilities = [branch["rank"] <= 1 for (b, branch) in sol["branch"]]
        branch_feasbility = all(branch_feasibilities)
        bus_feasibilities = [bus["rank"] <= 1 for (b, bus) in sol["bus"]]
        bus_feasibility = all(bus_feasibilities)

        ac_feasibility = bus_feasibility && branch_feasbility
        sol["is_ac_feasible_if_radial"] = ac_feasibility
    end
end


""
function add_rank_voltage_variable!(sol, pm::AbstractUBFModels; tol = 1e-6)
    for (nw, network) in pm._PMs.ref[:nw]
        buses       = _PMs.ref(pm, nw, :bus)
        for (b, bus) in buses
            Wre, Wim = _make_hermitian_matrix_variable(sol["bus"]["$b"]["w"].values, sol["bus"]["$b"]["wr"].values, sol["bus"]["$b"]["wi"].values)
            W = Wre + im*Wim
            sol["bus"]["$b"]["rank"] = rank(W, tol)
        end
    end
end


""
function add_rank_current_variable!(sol, pm::AbstractUBFModels; tol = 1e-6)
    for (nw, network) in pm._PMs.ref[:nw]
        branches       = _PMs.ref(pm, nw, :branch)
        for (b, branch) in branches
            CCre, CCim = _make_hermitian_matrix_variable(sol["branch"]["$b"]["cc"].values, sol["branch"]["$b"]["ccr"].values, sol["branch"]["$b"]["cci"].values)
            CC = CCre + im*CCim
            sol["branch"]["$b"]["CC_rank"] = rank(CC, tol)
        end
    end
end


""
function add_rank_branch_flow!(sol, pm::AbstractUBFModels; tol = 1e-6)
    for (nw, network) in pm._PMs.ref[:nw]
        buses       = _PMs.ref(pm, nw, :bus)
        for (b, branch) in _PMs.ref(pm, nw, :branch)
            g_fr = branch["g_fr"].values
            b_fr = branch["b_fr"].values
            y_fr = g_fr + im* b_fr

            fbus = branch["f_bus"]
            Wre, Wim = _make_hermitian_matrix_variable(sol["bus"]["$fbus"]["w"].values, sol["bus"]["$fbus"]["wr"].values, sol["bus"]["$fbus"]["wi"].values)
            W = Wre + im*Wim

            CCre, CCim = _make_hermitian_matrix_variable(sol["branch"]["$b"]["cc"].values, sol["branch"]["$b"]["ccr"].values, sol["branch"]["$b"]["cci"].values)
            CC = CCre + im*CCim

            Pij = _make_full_matrix_variable(sol["branch"]["$b"]["pf"].values, sol["branch"]["$b"]["pf_lt"].values, sol["branch"]["$b"]["pf_ut"].values)
            Qij = _make_full_matrix_variable(sol["branch"]["$b"]["qf"].values, sol["branch"]["$b"]["qf_lt"].values, sol["branch"]["$b"]["qf_ut"].values)
            Sij = Pij + im*Qij
            Ssij = Sij - W'*y_fr'

            f = [W Ssij; Ssij' CC]
            sol["branch"]["$b"]["rank"] = rank(f, tol)
        end
    end
end


""
function add_original_variables!(sol, pm::_PMs.AbstractPowerModel)
    if !haskey(pm.setting, "output") || !haskey(pm.setting["output"], "branch_flows") || pm.setting["output"]["branch_flows"] == false
        Memento.error(_LOGGER, "deriving the original variables requires setting: branch_flows => true")
    end

    for (nw, network) in pm._PMs.ref[:nw]
        #find rank-1 starting points
        ref_buses   = _find_ref_buses(pm, nw)
        #TODO develop code to start with any rank-1 W variable
        buses       = _PMs.ref(pm, nw, :bus)
        arcs        = _PMs.ref(pm, nw, :arcs)
        branches    = _PMs.ref(pm, nw, :branch)
        #define sets to explore
        all_bus_ids             = Set([b for (b, bus)      in _PMs.ref(pm, nw, :bus)])
        all_arc_from_ids        = Set([(l,i,j) for (l,i,j) in _PMs.ref(pm, nw, :arcs_from)])
        all_arc_to_ids          = Set([(l,i,j) for (l,i,j) in _PMs.ref(pm, nw, :arcs_to)])
        all_branch_ids          = Set([l for (l,i,j)       in _PMs.ref(pm, nw, :arcs_from)])
        visited_arc_from_ids    = Set()
        visited_arc_to_ids      = Set()
        visited_bus_ids         = Set()
        visited_branch_ids      = Set()

        for b in ref_buses
            sol["bus"]["$b"]["va"] = [0, -2*pi/3, 2*pi/3] #TODO support arbitrary angles at the reference bus
            sol["bus"]["$b"]["vm"] = _PMs.ref(pm, nw, :bus, b)["vm"].values
            push!(visited_bus_ids, b)
        end

        tt = 0
        while visited_branch_ids != all_branch_ids && visited_bus_ids != all_bus_ids
            tt = tt+1
            if tt >10000
                break
            end

            remaining_arc_from_ids = setdiff(all_arc_from_ids, visited_arc_from_ids)
            remaining_arc_to_ids = setdiff(all_arc_to_ids, visited_arc_to_ids)

            candidate_arcs_from = [(l,i,j) for (l,i,j) in remaining_arc_from_ids if i in visited_bus_ids && !(j in visited_bus_ids)]
            candidate_arcs_to   = [(l,i,j) for (l,i,j) in remaining_arc_to_ids   if i in visited_bus_ids && !(j in visited_bus_ids)]

            if !isempty(candidate_arcs_from)
                (l,i,j) = arc = candidate_arcs_from[1]
                g_fr = branches[l]["g_fr"].values
                b_fr = branches[l]["b_fr"].values
                y_fr = g_fr + im* b_fr
                g_to = branches[l]["g_to"].values
                b_to = branches[l]["b_to"].values
                y_to = g_to + im* b_to
                r = branches[l]["br_r"].values
                x = branches[l]["br_x"].values
                z = (r + im*x)
                Ui = sol["bus"]["$i"]["vm"].*exp.(im*sol["bus"]["$i"]["va"])

                Pij = _make_full_matrix_variable(sol["branch"]["$l"]["pf"].values, sol["branch"]["$l"]["pf_lt"].values, sol["branch"]["$l"]["pf_ut"].values)
                Qij = _make_full_matrix_variable(sol["branch"]["$l"]["qf"].values, sol["branch"]["$l"]["qf_lt"].values, sol["branch"]["$l"]["qf_ut"].values)
                Sij = Pij + im*Qij

                Ssij = Sij - Ui*Ui'*y_fr'
                Isij = (1/tr(Ui*Ui'))*(Ssij')*Ui
                Uj = Ui - z*Isij
                Iij = Isij + y_fr*Ui

                Isji = -Isij
                Iji = Isji + y_to*Uj

                sol["bus"]["$j"]["vm"] = abs.(Uj)
                sol["bus"]["$j"]["va"] = _wrap_to_pi(angle.(Uj))

                sol["branch"]["$l"]["cfm"] = abs.(Iij)
                sol["branch"]["$l"]["cfa"] = _wrap_to_pi(angle.(Iij))
                sol["branch"]["$l"]["ctm"] = abs.(Iji)
                sol["branch"]["$l"]["cta"] = _wrap_to_pi(angle.(Iji))
                #
                push!(visited_arc_from_ids, arc)
                push!(visited_branch_ids, l)
                push!(visited_bus_ids, j)

            elseif !isempty(candidate_arcs_to)
                (l,i,j) = arc = candidate_arcs_to[1]
                g_fr = branches[l]["g_to"].values
                b_fr = branches[l]["b_to"].values
                y_fr = g_fr + im* b_fr
                g_to = branches[l]["g_fr"].values
                b_to = branches[l]["b_fr"].values
                y_to = g_to + im* b_to
                r = branches[l]["br_r"].values
                x = branches[l]["br_x"].values
                z = (r + im*x)
                Ui = sol["bus"]["$i"]["vm"].*exp.(im*sol["bus"]["$i"]["va"])

                Pij = _make_full_matrix_variable(sol["branch"]["$l"]["pt"].values, sol["branch"]["$l"]["pt_lt"].values, sol["branch"]["$l"]["pt_ut"].values)
                Qij = _make_full_matrix_variable(sol["branch"]["$l"]["qt"].values, sol["branch"]["$l"]["qt_lt"].values, sol["branch"]["$l"]["qt_ut"].values)
                Sij = Pij + im*Qij

                Ssij = Sij - Ui*Ui'*y_fr'
                Isij = (1/tr(Ui*Ui'))*(Ssij')*Ui
                Uj = Ui - z*Isij
                Iij = Isij + y_fr*Ui

                Isji = -Isij
                Iji = Isji + y_to*Uj


                sol["bus"]["$j"]["vm"] = abs.(Uj)
                sol["bus"]["$j"]["va"] = _wrap_to_pi(angle.(Uj))

                sol["branch"]["$l"]["ctm"] = abs.(Iij)
                sol["branch"]["$l"]["cta"] = _wrap_to_pi(angle.(Iij))
                sol["branch"]["$l"]["cfm"] = abs.(Iji)
                sol["branch"]["$l"]["cfa"] = _wrap_to_pi(angle.(Iji))
                #
                push!(visited_arc_to_ids, arc)
                push!(visited_branch_ids, l)
                push!(visited_bus_ids, j)

            else #in case you have loops or multiple reference buses
                candidate_arcs = [(l,i,j) for (l,i,j) in remaining_arc_from_ids if i in visited_bus_ids && j in visited_bus_ids]
                (l,i,j) = arc = candidate_arcs[1]
                Sij = sol["branch"]["$l"]["pf"] + im* sol["branch"]["$l"]["qf"]
                Sji = sol["branch"]["$l"]["pt"] + im* sol["branch"]["$l"]["qt"]
                Ui = sol["bus"]["$i"]["vm"].*exp.(im*sol["bus"]["$i"]["va"])
                Uj = sol["bus"]["$j"]["vm"].*exp.(im*sol["bus"]["$j"]["va"])

                Iij = Sij./Ui
                Iji = Sji./Uj
                sol["branch"]["$l"]["cfm"] = abs.(Iij)
                sol["branch"]["$l"]["cfa"] = _wrap_to_pi(angle.(Iij))
                sol["branch"]["$l"]["ctm"] = abs.(Iji)
                sol["branch"]["$l"]["cta"] = _wrap_to_pi(angle.(Iji))
                push!(visited_arc_from_ids, arc)
                push!(visited_arc_to_ids, (l,j,i))
                push!(visited_branch_ids, l)
            end
        end
    end
end


"solution builder for minimum load delta problem (load shed)"
function solution_mld!(pm::_PMs.AbstractPowerModel, sol::Dict{String,Any})
    _PMs.add_setpoint_bus_voltage!(sol, pm)
    _PMs.add_setpoint_generator_power!(sol, pm)
    _PMs.add_setpoint_storage!(sol, pm)
    _PMs.add_setpoint_branch_flow!(sol, pm)

    add_setpoint_bus_status!(sol, pm)
    add_setpoint_load!(sol, pm)
    add_setpoint_shunt!(sol, pm)
    add_setpoint_generator_status!(sol, pm)
    add_setpoint_storage_status!(sol, pm)
end


"solution builder for branch-flow minimum load delta problem (load shed)"
function solution_mld_bf!(pm::_PMs.AbstractPowerModel, sol::Dict{String,Any})
    add_setpoint_bus_voltage!(sol, pm)
    _PMs.add_setpoint_generator_power!(sol, pm)
    add_setpoint_branch_flow!(sol, pm)
    add_setpoint_branch_current!(sol, pm)
    _PMs.add_setpoint_dcline_flow!(sol, pm)

    _PMs.add_dual_kcl!(sol, pm)
    _PMs.add_dual_sm!(sol, pm) # Adds the duals of the transmission lines' thermal limits.

    add_setpoint_bus_status!(sol, pm)
    add_setpoint_load!(sol, pm)
    add_setpoint_shunt!(sol, pm)
    add_setpoint_generator_status!(sol, pm)
    add_setpoint_storage_status!(sol, pm)

    if haskey(pm.setting, "output") && haskey(pm.setting["output"], "original_variables") && pm.setting["output"]["original_variables"] == true
        add_rank!(sol, pm)
        add_is_ac_feasible!(sol, pm)
        add_original_variables!(sol, pm)
    end
end


"add load setpoints for load shed problem"
function add_setpoint_load!(sol, pm::_PMs.AbstractPowerModel)
    _PMs.add_setpoint!(sol, pm, "load", "pd", :z_demand; scale = (x,item,cnd) -> x*item["pd"], conductorless=true)
    _PMs.add_setpoint!(sol, pm, "load", "qd", :z_demand; scale = (x,item,cnd) -> x*item["qd"], conductorless=true)
    _PMs.add_setpoint!(sol, pm, "load", "status", :z_demand; default_value = (item) -> if (item["status"] == 0) 0.0 else 1.0 end, conductorless=true)
end


"add shunt setpoints for load shed problem"
function add_setpoint_shunt!(sol, pm::_PMs.AbstractPowerModel)
    _PMs.add_setpoint!(sol, pm, "shunt", "gs", :z_shunt; scale = (x,item,cnd) -> x.*item["gs"], conductorless=true)
    _PMs.add_setpoint!(sol, pm, "shunt", "bs", :z_shunt; scale = (x,item,cnd) -> x.*item["bs"], conductorless=true)
    _PMs.add_setpoint!(sol, pm, "shunt", "status", :z_shunt; default_value = (item) -> if (item["status"] == 0) 0.0 else 1.0 end, conductorless=true)
end


"add bus statuses for load shed problem"
function add_setpoint_bus_status!(sol, pm::_PMs.AbstractPowerModel)
   _PMs.add_setpoint!(sol, pm, "bus", "status", :z_voltage; status_name="bus_type", inactive_status_value=4, default_value = (item) -> if item["bus_type"] == 4 0.0 else 1.0 end, conductorless=true)
end


"add generator statuses for load shed problem"
function add_setpoint_generator_status!(sol, pm::_PMs.AbstractPowerModel)
    _PMs.add_setpoint!(sol, pm, "gen", "gen_status", :z_gen; status_name="gen_status", default_value = (item) -> if (item["gen_status"] == 0) 0.0 else 1.0 end, conductorless=true)
end


"add storage statuses for load shed problem"
function add_setpoint_storage_status!(sol, pm::_PMs.AbstractPowerModel)
    _PMs.add_setpoint!(sol, pm, "storage", "status", :z_storage; default_value = (item) -> if (item["status"] == 0) 0.0 else 1.0 end, conductorless=true)
end
