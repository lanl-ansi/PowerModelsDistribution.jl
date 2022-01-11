# issues to launch
# - incorrect setup of short circuit impedance -> neutral row is all zeros, therefore impedance is not invertible
# - voltage source (gen 2) has g=0, qg=0 --> yprim=0 which causes singularity, voltage source does not have pg,qg keys anyways
# test_gen_3ph_wye.dss --> branch 2 (voltage source internal branch) has br_r,br_x=0

using Pkg
Pkg.activate("./scripts")
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


# case_file = joinpath(pwd(), "test/data/opendss/case_mxshunt_2.dss")                   # not very exact
# case_file = joinpath(pwd(), "test/data/opendss/case_mxshunt.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case2_diag.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_basefreq.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_battery.dss")           # not very exact
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_cap.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_isc.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_prop-order.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_pv.dss")                # not very exact
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced_switch.dss")            # yields NAN values
# case_file = joinpath(pwd(), "test/data/opendss/case3_balanced.dss")                   # error: requires load profile
# case_file = joinpath(pwd(), "test/data/opendss/case3_delta_gens.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case3_lm_1230.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case3_lm_models.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case3_unbalanced_1phase-pv.dss")       # not very exact
# case_file = joinpath(pwd(), "test/data/opendss/case3_unbalanced_assym_swap.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case3_unbalanced_delta_loads.dss")     # not very exact
# case_file = joinpath(pwd(), "test/data/opendss/case3_unbalanced_missingedge.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case3_unbalanced_switch.dss")          # yields NAN values
# case_file = joinpath(pwd(), "test/data/opendss/case3_unbalanced.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case4_phase_drop.dss")
# case_file = joinpath(pwd(), "test/data/opendss/case5_phase_drop.dss")

# case_file = joinpath(pwd(), "test/data/opendss/IEEE13_CapControl.dss")            # error
# case_file = joinpath(pwd(), "test/data/opendss/IEEE13_RegControl.dss")            # not exact
# case_file = joinpath(pwd(), "test/data/opendss/IEEE13_test_controls.dss")         # not exact

# case_file = joinpath(pwd(), "test/data/opendss/loadparser_error.dss")             # error: I think errors are expected for each load in thus case
# case_file = joinpath(pwd(), "test/data/opendss/loadparser_warn_model.dss")        # opendss does not run
# case_file = joinpath(pwd(), "test/data/opendss/test2_master.dss");                # too many errors

# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_2w_dy_lag.dss")
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_2w_dy_lead_small_series_impedance.dss")
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_2w_dy_lead.dss")
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_2w_yy_bank.dss")          # error fixed: xmfr leadlag --> spread into transformers leadlag
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_2w_yy_oltc.dss")
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_2w_yy.dss")
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_3w_dyy_1.dss")
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_3w_dyy_2.dss")            # error fixed: ODSS cannot handle xht=0 --> xht=0.1
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_3w_dyy_3_loadloss.dss")   # error fixed: xhl=0 xht=0 xlt=0 --> xhl=0.1 xht=0.1 xlt=0.1 
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_3w_dyy_3.dss")            # error fixed: xhl=0 xht=0 xlt=0 --> xhl=0.1 xht=0.1 xlt=0.1 
# case_file = joinpath(pwd(), "test/data/opendss/ut_trans_3w_dyy_basetest.dss")     # error fixed: xhl=0 xht=0 xlt=0 --> xhl=0.1 xht=0.1 xlt=0.1 
# case_file = joinpath(pwd(), "test/data/opendss/virtual_sourcebus.dss")



run_dss(open(f->read(f, String), case_file))

##
eng =_PMD.parse_file(case_file)

eng["voltage_source"]["source"]["rs"][4,4] = eng["voltage_source"]["source"]["rs"][1,1]
eng["voltage_source"]["source"]["rs"][1:3,4] .= eng["voltage_source"]["source"]["rs"][1,2]
eng["voltage_source"]["source"]["rs"][4,1:3] .= eng["voltage_source"]["source"]["rs"][1,2]
eng["voltage_source"]["source"]["xs"][4,4] = eng["voltage_source"]["source"]["xs"][1,1]
eng["voltage_source"]["source"]["xs"][1:3,4] .= eng["voltage_source"]["source"]["xs"][1,2]
eng["voltage_source"]["source"]["xs"][4,1:3] .= eng["voltage_source"]["source"]["xs"][1,2]


math = transform_data_model(eng;kron_reduce=false)

_PMD.add_start_vrvi!(math)
v_start = _PMD._bts_to_start_voltage(math)

result = compute_pf(math)

# @show result["iterations"]

# @show result["Yf"]

# basekv = math["settings"]["vbases_default"]["2"] * 1E3
# baseMVA = math["settings"]["sbase_default"] * 1E6
# Zbase = (400/sqrt(3))^2 / 1000
# (result["Yv"].L * result["Yv"].U) * Zbase


# G = [166807968.9    15223203.23    15223203.23   -166807968.9   -15223203.23   -15223203.23 
#     15223203.23    166807968.9    15223203.23   -15223203.23   -166807968.9   -15223203.23 
#     15223203.23    15223203.23    166807968.9   -15223203.23   -15223203.23   -166807968.9 
#    -166807968.9   -15223203.23   -15223203.23    166807968.9    15223203.23    15223203.23 
#    -15223203.23   -166807968.9   -15223203.23    15223203.23    166807968.9    15223203.23 
#    -15223203.23   -15223203.23   -166807968.9    15223203.23    15223203.23    166807968.9 ]

# B = [-601480417.1     4858645.53     4858645.53    601480417.1    -4858645.53    -4858645.53 
#     4858645.53   -601480417.1     4858645.53    -4858645.53    601480417.1    -4858645.53 
#     4858645.53     4858645.53   -601480417.1    -4858645.53    -4858645.53    601480417.1 
#     601480417.1    -4858645.53    -4858645.53   -601480417.1     4858645.53     4858645.53 
#     -4858645.53    601480417.1    -4858645.53     4858645.53   -601480417.1     4858645.53 
#     -4858645.53    -4858645.53    601480417.1     4858645.53     4858645.53   -601480417.1 ]

# Yprim_odss = G + B*im


# Yprim_pu = [
#     8.88799e9-3.20589e10im   8.03468e8+2.79172e8im    8.03468e8+2.79172e8im    8.03468e8+2.79172e8im   -8.88799e9+3.20589e10im  -8.03468e8-2.79172e8im   -8.03468e8-2.79172e8im   -8.03468e8-2.79172e8im
#     8.03468e8+2.79172e8im    8.88799e9-3.20589e10im   8.03468e8+2.79172e8im    8.03468e8+2.79172e8im   -8.03468e8-2.79172e8im   -8.88799e9+3.20589e10im  -8.03468e8-2.79172e8im   -8.03468e8-2.79172e8im
#     8.03468e8+2.79172e8im    8.03468e8+2.79172e8im    8.88799e9-3.20589e10im   8.03468e8+2.79172e8im   -8.03468e8-2.79172e8im   -8.03468e8-2.79172e8im   -8.88799e9+3.20589e10im  -8.03468e8-2.79172e8im
#     8.03468e8+2.79172e8im    8.03468e8+2.79172e8im    8.03468e8+2.79172e8im    8.88799e9-3.20589e10im  -8.03468e8-2.79172e8im   -8.03468e8-2.79172e8im   -8.03468e8-2.79172e8im   -8.88799e9+3.20589e10im
#    -8.88799e9+3.20589e10im  -8.03468e8-2.79172e8im   -8.03468e8-2.79172e8im   -8.03468e8-2.79172e8im    8.88799e9-3.20589e10im   8.03468e8+2.79172e8im    8.03468e8+2.79172e8im    8.03468e8+2.79172e8im
#    -8.03468e8-2.79172e8im   -8.88799e9+3.20589e10im  -8.03468e8-2.79172e8im   -8.03468e8-2.79172e8im    8.03468e8+2.79172e8im    8.88799e9-3.20589e10im   8.03468e8+2.79172e8im    8.03468e8+2.79172e8im
#    -8.03468e8-2.79172e8im   -8.03468e8-2.79172e8im   -8.88799e9+3.20589e10im  -8.03468e8-2.79172e8im    8.03468e8+2.79172e8im    8.03468e8+2.79172e8im    8.88799e9-3.20589e10im   8.03468e8+2.79172e8im
#    -8.03468e8-2.79172e8im   -8.03468e8-2.79172e8im   -8.03468e8-2.79172e8im   -8.88799e9+3.20589e10im   8.03468e8+2.79172e8im    8.03468e8+2.79172e8im    8.03468e8+2.79172e8im    8.88799e9-3.20589e10im]

# Yprim = Yprim_pu / Zbase
# phases_inds = [collect(1:3);collect(5:7)]
# Yprim_abc = Yprim[phases_inds, phases_inds]

# Z_odss = inv(Yprim_odss[1:3,1:3])
# Zprim = inv(Yprim_abc[1:3,1:3])
# real(Z_odss) - real(Zprim)
# imag(Z_odss) - imag(Zprim)


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
# va = [(i, bus_data["va"]) for (i,bus_data) in result["solution"]["bus"]]


