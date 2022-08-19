using Pkg
using PowerModelsDistribution
using JSON

# Pkg.activate("./examples")

const PMD = PowerModelsDistribution

include("Native Power Flow - Data Cleaning Functions.jl")
case_path = "examples/native_pf_testcases"
solution_dir = "examples/native_pf_testcases/solutions"

##

# ## 3 wire   ---  get these files from  https://github.com/sanderclaeys/DistributionTestCases.jl
case = "ieee13_pmd"
# case = "ieee34_pmd"
# case = "ieee123_pmd"

data_eng = parse_file("$case_path/$case.dss", transformations=[transform_loops!])
data_eng["is_kron_reduced"] = true
data_eng["settings"]["sbase_default"] = 1
vsource_correction_to_3w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false, phase_project=false)
sourcebus_voltage_vector_correction_3wire!(data_math, explicit_neutral=false)
update_math_model_3wire!(data_math)
res = PowerModelsDistribution.compute_pf(data_math; explicit_neutral=false)

sol_dss = open("$solution_dir/$case.json", "r") do f
    JSON.parse(f)
end
sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)


## 3wire and 4wire native_pf_testcases

case = "test_trans_dy"
data_eng = parse_file("$case_path/$case.dss", transformations=[transform_loops!])
vsource_correction_to_4w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false)
res = PowerModelsDistribution.compute_pf(data_math; explicit_neutral=true)

sol_dss = open("$solution_dir/$case.json", "r") do f
    JSON.parse(f)
end
sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)


case = "test_trans_yy"
data_eng = parse_file("$case_path/$case.dss", transformations=[transform_loops!])
vsource_correction_to_4w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false)
res = PowerModelsDistribution.compute_pf(data_math; explicit_neutral=true)

sol_dss = open("$solution_dir/$case.json", "r") do f
    JSON.parse(f)
end
sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)



case = "test_trans_dy_3w"
data_eng = parse_file("$case_path/$case.dss", transformations=[transform_loops!])
data_eng["is_kron_reduced"] = true
vsource_correction_to_3w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false)
sourcebus_voltage_vector_correction_3wire!(data_math, explicit_neutral=false)
update_math_model_3wire!(data_math)
res = PowerModelsDistribution.compute_pf(data_math; explicit_neutral=false)

sol_dss = open("$solution_dir/$case.json", "r") do f
    JSON.parse(f)
end
sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)



case = "test_trans_yy_3w"
data_eng = parse_file("$case_path/$case.dss", transformations=[transform_loops!])
data_eng["is_kron_reduced"] = true
vsource_correction_to_3w!(data_eng)
data_math = transform_data_model(data_eng;kron_reduce=false)
sourcebus_voltage_vector_correction_3wire!(data_math, explicit_neutral=false)
update_math_model_3wire!(data_math)
res = PowerModelsDistribution.compute_pf(data_math; explicit_neutral=false)

sol_dss = open("$solution_dir/$case.json", "r") do f
    JSON.parse(f)
end
sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)


# obtain solution from dss


##
# run the unit tests with KLU, and with original factorize
# accessing OpenDSSDirect for results consistency
# validating the results are consistent with OpenDSS
# double check what is actually happeing with the aux current and voltage variables in ideal transformers
# fixing the unit tests
# double checking KLU settings (default or not?) with respect to OpenDSS settings


## solution builder - > make sure load currents are in dict (phase->cr, ci) format

##
# data_dir = "test/data/en_validation_case_data"
# cases = [x[1:end-4] for x in readdir(data_dir) if endswith(x, ".dss")]

# case = "test_load_1ph_delta_cp"
# case_path = "$data_dir/$case.dss"

# data_eng = parse_file(case_path, transformations=[transform_loops!])
# vsource_correction_to_4w!(data_eng)

# data_math = transform_data_model(data_eng;kron_reduce=false)

# res = compute_pf(data_math; explicit_neutral=true)

# ##

# vm = [v for (i,v) in sort(res["solution"]["bus"]["1"]["vm"])]  # load bus
# va = [v for (i,v) in sort(res["solution"]["bus"]["1"]["va"])]
# v = vm .* exp.(va*im)

# # c2r = res["solution"]["load"]["2"]["cdr"]
# # c2i = res["solution"]["load"]["2"]["cdi"]
# # c2 = c2r .+ im*c2i
# # S2 = (v[2] - v[3]) * c2[1]'

# c1r = res["solution"]["load"]["1"]["cdr"]
# c1i = res["solution"]["load"]["1"]["cdi"]
# c1 = c1r .+ im*c1i
# S1 = (v[1] - v[4]) * c1[1]'


# ##
# vm = [v for (i,v) in sort(res["solution"]["bus"]["3"]["vm"])]  # source bus
# va = [v for (i,v) in sort(res["solution"]["bus"]["3"]["va"])]
# v = vm .* exp.(va*im)

# cr_branch = res["solution"]["branch"]["2"]["cr"][1:4]
# ci_branch = res["solution"]["branch"]["2"]["ci"][1:4]
# c_banch = cr_branch .+ im*ci_branch

# using LinearAlgebra
# Ssource = diag(v * c_banch')
