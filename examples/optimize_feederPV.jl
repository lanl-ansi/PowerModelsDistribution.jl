import CSV
import JuMP
import Ipopt
import Pkg

using DataFrames
using Xpress
using Revise
using Printf

using PowerModelsDistribution
const PMD = PowerModelsDistribution

network_file = "test/data/opendss/123bus/IEEE123Master.dss"

data = parse_file(network_file, data_model = DSS)

apply_voltage_bounds!(data, vm_lb=0.97, vm_ub=1.03)

increase_factor = 4
Vulnerabele = ["60", "76", "97", "96", "114", "151"]
for (k, solar) in data["solar"]
    solar["pg"] = solar["pg"]*increase_factor
    solar["pg_ub"] = solar["pg_ub"]*increase_factor
    solar["pg_lb"] = solar["pg_lb"]*increase_factor
    solar["qg_ub"] = solar["qg_ub"]*increase_factor
    solar["qg_lb"] = solar["qg_lb"]*increase_factor
    if solar["bus"] in Vulnerabele
        solar["vulnerable"] = 1
    else
        solar["vulnerable"] = 0
    end
end



pm_type = LinDist3FlowPowerModel

pm = instantiate_mc_model(data, pm_type, build_mc_opf_curt)

# xpress_solver = optimizer_with_attributes(() -> Xpress.Optimizer(DEFAULTALG=2, PRESOLVE=0, logfile = "output.log"))
ipopt_solver = optimizer_with_attributes(() -> Ipopt.Optimizer())

res = optimize_model!(pm, optimizer=ipopt_solver)

sbase = data["settings"]["sbase_default"]
println('#', " \t", "Curtailed", " \t\t", "Available")
for (k, gen) in res["solution"]["gen"]
    pg = gen["pg"]*sbase
    pg_tot = sum(pg)
    
    pgvar = var(pm, :pg, parse(Int, k))
    if k == "17"
        Pglimit = 0
        continue
    else
        Pglimit = sum(JuMP.upper_bound(pgvar[i]) for i in eachindex(pgvar))*sbase
    end

    @printf("%s \t %.3f \t %.3f\n", k, Pglimit - pg_tot,  Pglimit)
end

Vub_ = 0
Vlb_ = Inf
for (k, bus) in res["solution"]["bus"]
    # vl = sqrt(minimum(bus["w"])) 
    # vh = sqrt(maximum(bus["w"]))
    if sqrt(minimum(bus["w"])) < Vlb_
        global Vlb_ = sqrt(minimum(bus["w"]))
    end
    if sqrt(maximum(bus["w"])) > Vub_
        global Vub_ = sqrt(maximum(bus["w"]))
    end
end
println("Maximum Voltage:  ", Vub_)
println("Minimum Voltage:  ", Vlb_)