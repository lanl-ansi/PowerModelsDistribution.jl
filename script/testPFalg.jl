using Pkg
projectFolder="C:\\CSIRO_Docs\\55_PMD_Fork\\PowerModelsDistribution.jl\\"
# issues to launch
# - incorrect setup of short circuit impedance -> neutral row is all zeros, therefore impedance is not invertible
# data_foler=projectFolder*"test\\data\\opendss\\"
# data_file_name="case3_unbalanced.dss"

data_foler=projectFolder*"test\\data\\en_validation_case_data\\"
data_file_name="test_gen_3ph_wye.dss"

cd(projectFolder)
Pkg.activate(".")
using PowerModelsDistribution
const _PMD=PowerModelsDistribution

I = [1 0 0 0; 0 1 0 0 ; 0 0 1 0; 0 0 0 1]

eng=_PMD.parse_file(data_foler*data_file_name)
math = transform_data_model(eng;kron_reduce=false)

math["branch"]["2"]["br_r"] = 1E-3.*I
math["branch"]["2"]["br_x"] = 1E-3.*I

_PMD.add_start_vrvi!(math)
v_start = _PMD._bts_to_start_voltage(math)

# result=compute_pf(math)

pfd = PowerFlowData(math, v_start)
stat_tol=1E-8
(Uv, status, its, stat) = _PMD._compute_Uv(pfd, max_iter=20, stat_tol=stat_tol)

sol= build_solution(pfd, Uv)



