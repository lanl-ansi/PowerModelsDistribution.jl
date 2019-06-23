"Adds sequence components of the voltage to the solution dict."
function get_solution_vm_all(pm::PMs.GenericPowerModel, sol::Dict{String,<:Any})
    PMs.get_solution(pm, sol)
    if !PMs.ismultinetwork(pm)
        get_solution_vm_all(pm, pm.cnw, sol)
    else
        for nw in PMs.nws(pm)
            get_solution_vm_all(pm, nw, sol[nw])
        end
    end
end

function get_solution_vm_all(pm::PMs.GenericPowerModel, nw::Int, sol_nw::Dict{String,<:Any})
    PMs.get_solution(pm, sol_nw)
    if !PMs.ismulticonductor(pm, nw=nw) && PMs.ref(pm, nw, :conductors)==3
        Memento.error(LOGGER, "Sequence components are only defined on three-phase networks.")
    end
    for bus_id in PMs.ids(pm, nw, :bus)
        vabc = [sol_nw["bus"]["$bus_id"]["vm"][c]*exp(im*sol_nw["bus"]["$bus_id"]["va"][c]) for c in 1:3]
        a = exp(im*2*pi/3)
        vpos = transpose(vabc)*[1, a, a^2]./3
        vneg = transpose(vabc)*[1, a^2, a]./3
        vzero = transpose(vabc)*[1, 1, 1]./3
        sol_nw["bus"]["$bus_id"]["vm_pos_seq"] = abs(vpos)
        sol_nw["bus"]["$bus_id"]["vm_neg_seq"] = abs(vneg)
        sol_nw["bus"]["$bus_id"]["vm_zero_seq"] = abs(vzero)
        sol_nw["bus"]["$bus_id"]["vuf"] = abs(vneg)/abs(vpos)
        sol_nw["bus"]["$bus_id"]["vm_ll"] = PMs.MultiConductorVector(abs.([1 -1 0; 0 1 -1; -1 0 1]*vabc)./sqrt(3))
    end
end
