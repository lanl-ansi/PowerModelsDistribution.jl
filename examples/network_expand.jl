import CSV
import JuMP
import Ipopt
import Pkg

using DataFrames
using Xpress
using Revise

using PowerModelsDistribution
const PMD = PowerModelsDistribution

# ENV["CPLEX_STUDIO_BINARIES"] = "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/"
# Pkg.add("CPLEX")
# Pkg.build("CPLEX")

using CPLEX

cplex_solver = optimizer_with_attributes(() -> CPLEX.Optimizer())
xpress_solver = optimizer_with_attributes(() -> Xpress.Optimizer(DEFAULTALG=2, PRESOLVE=0, logfile = "output.log"))


network_file = "test/data/opendss/IEEE13_CDPSM.dss";
expansion_file = "test/data/ne_json/expansion13.json"

data = parse_file(network_file)
expand_mods = parse_file(expansion_file)

PMD._IM.update_data!(data, expand_mods)

pm_type = NFAUPowerModel

pm = instantiate_mc_model(data, pm_type, build_mn_mc_mld_simple)
expand_cost = objective_ne(pm)
JuMP.@constraint(pm.model, expand_cost <= 1000)
objective_mc_min_fuel_cost(pm)

res = optimize_model!(pm, optimizer=cplex_solver)