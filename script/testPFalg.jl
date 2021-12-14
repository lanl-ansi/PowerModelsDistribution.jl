# issues to launch
# - incorrect setup of short circuit impedance -> neutral row is all zeros, therefore impedance is not invertible
# data_foler=projectFolder*"test\\data\\opendss\\"
# data_file_name="case3_unbalanced.dss"
# - voltage source (gen 2) has g=0, qg=0 --> yprim=0 which causes singularity, voltage source does not have pg,qg keys anyways
# test_gen_3ph_wye.dss --> branch 2 (voltage source internal branch) has br_r,br_x=0

using Pkg
Pkg.activate("./scripts")
# Pkg.develop("./")
using PowerModelsDistribution
const _PMD=PowerModelsDistribution


# case_file = joinpath(dirname(pathof(PowerModelsDistribution)), "../test/data/en_validation_case_data/test_gen_3ph_wye.dss")
case_file = "/Users/hei06j/Documents/repositories/remote/PowerModelsDistribution.jl/test/data/en_validation_case_data/test_trans_dy.dss"

##
I = [1 0 0 0; 0 1 0 0 ; 0 0 1 0; 0 0 0 1]
# I = [1 0 0; 0 1 0 ; 0 0 1]

eng =_PMD.parse_file(case_file)
math = transform_data_model(eng;kron_reduce=false)

# math["branch"]["2"]["br_r"] = 1E-3.*I
# math["branch"]["2"]["br_x"] = 1E-3.*I

math["branch"]["6"]["br_r"] = 1E-6.*I
math["branch"]["6"]["br_x"] = 1E-6.*I

_PMD.add_start_vrvi!(math)
v_start = _PMD._bts_to_start_voltage(math)

result=compute_pf(math)

# pfd = PowerFlowData(math, v_start)
# stat_tol=1E-8
# (Uv, status, its, stat) = _PMD._compute_Uv(pfd, max_iter=20, stat_tol=stat_tol)

# sol= build_solution(pfd, Uv)



