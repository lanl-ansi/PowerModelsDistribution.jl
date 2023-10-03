import CSV
import JuMP
import Ipopt
import Pkg

using DataFrames
using Revise
using Printf
using Juniper
using SCIP

using PowerModelsDistribution
const PMD = PowerModelsDistribution

network_file = "/Users/azamzam/dev/Julia/power-water-expansion/PowerModelsDistribution.jl/test/data/opendss/case3_unbalanced.dss"
expansion_file = "/Users/azamzam/dev/Julia/power-water-expansion/PowerModelsDistribution.jl/test/data/ne_json/expansion3u.json"

data = parse_file(network_file, data_model = DSS)
expand_mods = parse_file(expansion_file)
PMD._IM.update_data!(data, expand_mods)

apply_voltage_bounds!(data, vm_lb=0.97, vm_ub=1.03)

data_mn = make_multinetwork(data)

interval_len = Dict("interval_len"=> 24)
PMD._IM.update_data!(data_mn, interval_len)

pm_type = ACPUPowerModel

data_mn_math = transform_data_model(data_mn; kron_reduce=true)

pm = instantiate_mc_model(data_mn_math, pm_type, build_mn_mc_mld_multi_scenario_ne)
expand_cost = objective_ne(pm)
JuMP.@constraint(pm.model, expand_cost <= 100000)
objective_mc_min_fuel_cost(pm)

ipopt_solver = optimizer_with_attributes(() -> Ipopt.Optimizer())
juniper_optimizer = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>ipopt_solver, "time_limit"=>600)
scip_optimizer = optimizer_with_attributes(SCIP.Optimizer, "limits/gap"=>0.05)

res = optimize_model!(pm, optimizer=juniper_optimizer)
