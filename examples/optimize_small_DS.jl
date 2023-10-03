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

network_file = "test/data/opendss/case3_unbalanced.dss"

data = parse_file(network_file, data_model = DSS)

increase_factor = 4

pm_type = ACPUPowerModel

data_math = transform_data_model(data; kron_reduce=true)

apply_voltage_bounds!(data_math, vm_lb=0.97, vm_ub=1.03)

pm = instantiate_mc_model(data_math, pm_type, build_mc_opf)


# xpress_solver = optimizer_with_attributes(() -> Xpress.Optimizer(DEFAULTALG=2, PRESOLVE=0, logfile = "output.log"))
ipopt_solver = optimizer_with_attributes(() -> Ipopt.Optimizer())

res = optimize_model!(pm, optimizer=ipopt_solver)
