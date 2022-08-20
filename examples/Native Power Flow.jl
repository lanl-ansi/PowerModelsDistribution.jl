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
# accessing OpenDSSDirect for results consistency
# validating the results are consistent with OpenDSS
# double check what is actually happeing with the aux current and voltage variables in ideal transformers
# fixing the unit tests
# double checking KLU settings (default or not?) with respect to OpenDSS settings


## solution builder - > make sure load currents are in dict (phase->cr, ci) format