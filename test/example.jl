using Revise

using PowerModelsDistribution
const PMD = PowerModelsDistribution

import InfrastructureModels

import JuMP
import Ipopt

import JSON


using LinearAlgebra

pmd_path = joinpath(dirname(pathof(PowerModelsDistribution)), ".")

ipopt_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
ipopt_solver_adaptive = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "mu_strategy"=>"adaptive")

case3_balanced = parse_file("./test/data/opendss/case3_balanced_nobattery.dss")
modifications = parse_file("./test/data/ne_json/case3_balanced_ne.json")
modifications["storage_ne"]["s_ne1"]["configuration"] = PMD.WYE
modifications["storage_ne"]["s_ne1"]["status"] = PMD.ENABLED

PMD._IM.update_data!(case3_balanced, modifications)

# case3_balanced_battery = parse_file("./test/data/opendss/case3_balanced_battery.dss")

data_mn_nobatt = make_multinetwork(case3_balanced)

# data_mn_batt   = make_multinetwork(case3_balanced_battery)

res1 = solve_mn_mc_mld_simple_ne(data_mn_nobatt, NFAUPowerModel, ipopt_solver, relax_integrality=true)
# res2 = solve_mn_mc_mld_simple(data_mn_batt, NFAUPowerModel, ipopt_solver)


