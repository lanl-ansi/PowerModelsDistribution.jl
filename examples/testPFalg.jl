using Pkg
Pkg.activate("./examples")
# Pkg.develop("./")
using PowerModelsDistribution
using OpenDSSDirect
const _PMD=PowerModelsDistribution
const _ODSS=OpenDSSDirect

function run_dss(filename)
    dss(filename)
    dss("""
    Set Toler=0.000000001
    // Dump Line.*  debug
    Show Voltages LN Nodes
    // Show Yprim
    // Show Y
    // Dump Debug
    """)
end

##

path = joinpath(dirname(pathof(PowerModelsDistribution)), "..")
cd(path)


### fix the issue with these cases
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_gen_1ph_wye.dss")
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_gen_1ph_delta.dss")
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_gen_3ph_wye.dss")
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_gen_3ph_delta.dss")
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_grounding.dss")
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_load_1ph_delta_cp.dss")
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_load_1ph_wye_cp.dss")
case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_load_3ph_delta_ci.dss")
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_load_3ph_delta_cp.dss")
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_load_3ph_delta_cz.dss")   # not very exact
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_load_3ph_wye_ci_balanced_load.dss")     # not very exact
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_load_3ph_wye_cp_balanced_load.dss")     # not very exact
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_load_3ph_wye_cz_balanced_load.dss")     # not very exact
# case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_trans_dy.dss")


# run_dss(open(f->read(f, String), case_file))

##
eng =_PMD.parse_file(case_file, transformations=[transform_loops!])

eng["voltage_source"]["source"]["rs"][4,4] = eng["voltage_source"]["source"]["rs"][1,1]
eng["voltage_source"]["source"]["rs"][1:3,4] .= eng["voltage_source"]["source"]["rs"][1,2]
eng["voltage_source"]["source"]["rs"][4,1:3] .= eng["voltage_source"]["source"]["rs"][1,2]
eng["voltage_source"]["source"]["xs"][4,4] = eng["voltage_source"]["source"]["xs"][1,1]
eng["voltage_source"]["source"]["xs"][1:3,4] .= eng["voltage_source"]["source"]["xs"][1,2]
eng["voltage_source"]["source"]["xs"][4,1:3] .= eng["voltage_source"]["source"]["xs"][1,2]


math = transform_data_model(eng;kron_reduce=false)
result = compute_pf(math)

vm = [(i, bus_data["vm"]) for (i,bus_data) in result["solution"]["bus"]]

va = Dict()
for (i,bus_data) in result["solution"]["bus"]
    va[i] = Dict()
    va[i]["va"] = Dict()
    for (phase, va_phase) in bus_data["va"]
        va[i]["va"][phase] = va_phase * 180/pi
    end
end
va


_PMD.solution_make_si(result["solution"], math)
