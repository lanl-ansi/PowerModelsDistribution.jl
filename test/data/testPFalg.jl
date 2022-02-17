using Pkg
# Pkg.activate("./")

Pkg.activate("./test")
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
# path = joinpath(dirname(pathof(PowerModelsDistribution)), "..")
path = pwd()
cd(path)

case_file = joinpath(pwd(), "test/data/en_validation_case_data/test_load_3ph_wye_exp.dss")  


### fix the issue with these cases          transformations=[transform_loops!]      causes inexactness
# case_file = joinpath(pwd(), "test/data/opendss/case_mxshunt_2.dss")                   # not exact    
# case_file = joinpath(pwd(), "test/data/opendss/case_mxshunt.dss")                     # ✓
# case_file = joinpath(pwd(), "test/data/opendss/case2_diag.dss")                       # ✓
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_basefreq.dss")          # ✓
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_battery.dss")           # not exact
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_cap.dss")               # ✓
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_isc.dss")               # ✓
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_prop-order.dss")        # ✓
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_pv.dss")                # not exact
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_switch.dss")            # ✓ slight difference: switch br_r and br_x are not parsed and used
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced.dss")                   # ✓ without load_profile (ERRORED: load_profile.csv path updated)
# case_file = joinpath(pwd(), "test/data/opendss/case3_delta_gens.dss")                 # ✓
# case_file = joinpath(pwd(), "test/data/opendss/case3_lm_1230.dss")                    # ✓  un-interesting?
# case_file = joinpath(pwd(), "test/data/opendss/case3_lm_models.dss")                  # ✓  dss file updated
# case_file = joinpath(pwd(), "test/data/opendss/case3_unbalanced_1phase-pv.dss")       # ✓
# case_file = joinpath(pwd(), "test/data/opendss/case3_unbalanced_assym_swap.dss")      # ✓
# case_file = joinpath(pwd(), "test/data/opendss/case3_unbalanced_delta_loads.dss")     # throw error in unittest: not exact:  each arm of delta load is a different load model: compute_pf does not support this
# case_file = joinpath(pwd(), "test/data/opendss/case3_unbalanced_missingedge.dss")     # ✓
# case_file = joinpath(pwd(), "test/data/opendss/case3_unbalanced_switch.dss")          # ✓ slight difference:  switches with same length as lines, but impedance is not considered
# case_file = joinpath(pwd(), "test/data/opendss/case3_unbalanced.dss")                 # ✓
# case_file = joinpath(pwd(), "test/data/opendss/case4_phase_drop.dss")                 # ✓
# case_file = joinpath(pwd(), "test/data/opendss/case5_phase_drop.dss")                 # ✓

# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_2w_dy_lag.dss")               # ✓ slight difference
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_2w_dy_lead_small_series_impedance.dss")  # ✓
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_2w_dy_lead.dss")              # ✓ slight difference
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_2w_yy_oltc.dss")              # ✓
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_2w_yy.dss")                   # ✓
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_3w_dyy_1.dss")                # ✓ slight difference

# case_file = joinpath(pwd(), "test/data/opendss/IEEE13_CapControl.dss")                # ? error in opendss
# case_file = joinpath(pwd(), "test/data/opendss/IEEE13_RegControl.dss")                # ? not exact
# case_file = joinpath(pwd(), "test/data/opendss/IEEE13_test_controls.dss")             # ? not exact
# case_file = joinpath(pwd(), "test/data/opendss/loadparser_error.dss")                 # ? error: I think errors are expected for each load in this case
# case_file = joinpath(pwd(), "test/data/opendss/loadparser_warn_model.dss")            # ? error in opendss
# case_file = joinpath(pwd(), "test/data/opendss/test2_master.dss");                    # ? too many errors
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_3w_dyy_2.dss")                # ✓ slight difference: (ERRORED: Zero Reactance specified for Transformer.tx1 --> fixed: ODSS cannot handle xht=0 : has to stay the same due to other unit tests
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_3w_dyy_3_loadloss.dss")       # ✓ slight difference: (ERRORED: Zero Reactance specified for Transformer.tx1 : has to stay the same due to other unit tests
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_3w_dyy_3.dss")                # ✓ slight difference: (ERRORED: Zero Reactance specified for Transformer.tx1 : has to stay the same due to other unit tests
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_3w_dyy_basetest.dss")         # ✓ slight difference: (ERRORED: Zero Reactance specified for Transformer.tx1 : has to stay the same due to other unit tests
# case_file = joinpath(pwd(), "test/data/opendss/virtual_sourcebus.dss")                # ✓ slight difference


run_dss(open(f->read(f, String), case_file))

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

# vm = Dict()
# for (i,bus_data) in result["solution"]["bus"]
#     vm[i] = Dict()
#     vm[i]["vm"] = Dict()
#     for (phase, vm_phase) in bus_data["vm"]
#         vm[i]["vm"][phase] = vm_phase * 400/sqrt(3)
#     end
# end
# vm


# va = Dict()
# for (i,bus_data) in result["solution"]["bus"]
#     va[i] = Dict()
#     va[i]["va"] = Dict()
#     for (phase, va_phase) in bus_data["va"]
#         va[i]["va"][phase] = va_phase * 180/pi
#     end
# end
# va


# _PMD.solution_make_si(result["solution"], math)
