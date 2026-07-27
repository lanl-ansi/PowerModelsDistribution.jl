# 2-bus
case2_diag = parse_file("../test/data/opendss/case2_diag.dss")
case2_mxshunt = parse_file("../test/data/opendss/case_mxshunt.dss")
case2_mxshunt_2 = parse_file("../test/data/opendss/case_mxshunt_2.dss")
case2_virtual_sourcebus = parse_file("../test/data/opendss/virtual_sourcebus.dss")

# 3-bus balanced
case3_balanced = parse_file("../test/data/opendss/case3_balanced.dss")
case3_balanced_cap = parse_file("../test/data/opendss/case3_balanced_cap.dss")
case3_balanced_isc = parse_file("../test/data/opendss/case3_balanced_isc.dss")
case3_balanced_pv = parse_file("../test/data/opendss/case3_balanced_pv.dss")
case3_unbalanced_1phase_pv = parse_file("../test/data/opendss/case3_unbalanced_1phase-pv.dss")
case3_balanced_basefreq = parse_file("../test/data/opendss/case3_balanced_basefreq.dss")
case3_balanced_battery = parse_file("../test/data/opendss/case3_balanced_battery.dss")
case3_balanced_switch = parse_file("../test/data/opendss/case3_balanced_switch.dss")

# 3-bus unbalanced
case3_unbalanced = parse_file("../test/data/opendss/case3_unbalanced.dss")
case3_unbalanced_assym_swap = parse_file("../test/data/opendss/case3_unbalanced_assym_swap.dss")
case3_unbalanced_delta_loads = parse_file("../test/data/opendss/case3_unbalanced_delta_loads.dss")
case3_unbalanced_missingedge = parse_file("../test/data/opendss/case3_unbalanced_missingedge.dss")
case3_unbalanced_switch = parse_file("../test/data/opendss/case3_unbalanced_switch.dss")

# 3-bus load models
case3_lm_1230 = parse_file("../test/data/opendss/case3_lm_1230.dss")
case3_lm_models = parse_file("../test/data/opendss/case3_lm_models.dss")
case3_unbalanced_ZIPloads = parse_file("../test/data/opendss/case3_unbalanced_ZIPloads.dss")

# 3-bus delta generators
case3_delta_gens = parse_file("../test/data/opendss/case3_delta_gens.dss")

# 4-bus
case4_phase_drop = parse_file("../test/data/opendss/case4_phase_drop.dss")

# 5-bus
case5_phase_drop = parse_file("../test/data/opendss/case5_phase_drop.dss")

# 4-bus 2w transformer
ut_trans_2w_yy = parse_file("../test/data/opendss/ut_trans_2w_yy.dss")
ut_trans_2w_dy_lead = parse_file("../test/data/opendss/ut_trans_2w_dy_lead.dss")
ut_trans_2w_dy_lag = parse_file("../test/data/opendss/ut_trans_2w_dy_lag.dss")
ut_trans_2w_dy_lead_small_series_impedance = parse_file("../test/data/opendss/ut_trans_2w_dy_lead_small_series_impedance.dss", data_model=MATHEMATICAL)
ut_trans_2w_yy_bank = parse_file("../test/data/opendss/ut_trans_2w_yy_bank.dss")
ut_trans_2w_yy_unbanked = parse_file("../test/data/opendss/ut_trans_2w_yy_bank.dss"; bank_transformers=false)
ut_trans_2w_yy_oltc = parse_file("../test/data/opendss/ut_trans_2w_yy_oltc.dss")

# 4-bus 3w transformer
ut_trans_3w_dyy_1 = parse_file("../test/data/opendss/ut_trans_3w_dyy_1.dss")
ut_trans_3w_dyy_2 = parse_file("../test/data/opendss/ut_trans_3w_dyy_2.dss")
ut_trans_3w_dyy_3 = parse_file("../test/data/opendss/ut_trans_3w_dyy_3.dss")
ut_trans_3w_dyy_3_loadloss = parse_file("../test/data/opendss/ut_trans_3w_dyy_3_loadloss.dss")
trans_3w_center_tap = parse_file("../test/data/opendss/trans_3w_center_tap.dss")

# IEEE13
IEEE13_Assets = parse_file("../test/data/opendss/IEEE13_Assets.dss")
IEEE13_RegControl = parse_file("../test/data/opendss/IEEE13_RegControl.dss"; transformations=[remove_line_limits!, remove_transformer_limits!])
IEEE13_CapControl = parse_file("../test/data/opendss/IEEE13_CapControl.dss"; transformations=[remove_line_limits!, remove_transformer_limits!])

# explicit neutral
test_gen_3ph_wye = parse_file("../test/data/en_validation_case_data/test_gen_3ph_wye.dss", transformations=[remove_all_bounds!])
test_switch = parse_file("../test/data/en_validation_case_data/test_switch.dss", transformations=[remove_all_bounds!])
test_gen_1ph_wye = parse_file("../test/data/en_validation_case_data/test_gen_1ph_wye.dss", transformations=[remove_all_bounds!])
test_trans_dy = parse_file("../test/data/en_validation_case_data/test_trans_dy.dss", transformations=[remove_all_bounds!, transform_loops!])

# distribution transformer equivalent cases
dist_transformer = parse_file("../test/data/opendss/dist_transformer.dss")

# IEEE13 Feeder
IEEE13_Feeder_engr = parse_file("../test/data/opendss/ieee13_feeder.dss", multinetwork=true,
    time_series="daily", transformations=[remove_line_limits!, remove_transformer_limits!])

# RAVENS
ravens_case3_withgens = parse_file("../test/data/ravens/ravens_case3_withgens.json")
ravens_case3_withpvandstorage = parse_file("../test/data/ravens/ravens_case3_withpvandstorage.json")
ravens_case3_withsubxf = parse_file("../test/data/ravens/ravens_case3_withsubxf.json")
ravens_case3_withcap = parse_file("../test/data/ravens/ravens_case3_withcap.json")
ravens_test_switch_3w = parse_file("../test/data/ravens/ravens_test_switch_3w.json")
# ravens_test_switch_1w = parse_file("../test/data/ravens/ravens_test_switch_1w.json")
ravens_case3_withgens_mn = parse_file("../test/data/ravens/ravens_case3_withgens_mn.json")




# RAVENS-EXTRA
ravens_case2_diag = parse_file("../test/data/ravens/ravens_case2_diag.json")
ravens_case2_mxshunt = parse_file("../test/data/ravens/ravens_case_mxshunt.json")
ravens_case2_mxshunt_2 = parse_file("../test/data/ravens/ravens_case_mxshunt_2.json")
ravens_case2_virtual_sourcebus = parse_file("../test/data/ravens/ravens_virtual_sourcebus.json")

# 3-bus balanced
ravens_case3_balanced = parse_file("../test/data/ravens/ravens_case3_balanced.json")
ravens_case3_balanced = parse_file("../test/data/ravens/ravens_case3_balanced.json")
ravens_case3_balanced_cap = parse_file("../test/data/ravens/ravens_case3_balanced_cap.json")
ravens_case3_balanced_isc = parse_file("../test/data/ravens/ravens_case3_balanced_isc.json")
ravens_case3_balanced_pv = parse_file("../test/data/ravens/ravens_case3_balanced_pv.json")
ravens_case3_unbalanced_1phase_pv = parse_file("../test/data/ravens/ravens_case3_unbalanced_1phase-pv.json")
ravens_case3_balanced_basefreq = parse_file("../test/data/ravens/ravens_case3_balanced_basefreq.json")
ravens_case3_balanced_battery = parse_file("../test/data/ravens/ravens_case3_balanced_battery.json")
ravens_case3_balanced_switch = parse_file("../test/data/ravens/ravens_case3_balanced_switch.json")

# 3-bus unbalanced
ravens_case3_unbalanced = parse_file("../test/data/ravens/ravens_case3_unbalanced.json")
ravens_case3_unbalanced_assym_swap = parse_file("../test/data/ravens/ravens_case3_unbalanced_assym_swap.json")

ravens_case3_unbalanced_delta_loads = parse_file("../test/data/ravens/ravens_case3_unbalanced_delta_loads.json")
ravens_case3_unbalanced_delta_loads = parse_file("../test/data/ravens/ravens_case3_unbalanced_delta_loads.json")

ravens_case3_unbalanced_missingedge = parse_file("../test/data/ravens/ravens_case3_unbalanced_missingedge.json")
ravens_case3_unbalanced_switch = parse_file("../test/data/ravens/ravens_case3_unbalanced_switch.json")

# 3-bus load models
ravens_case3_lm_1230 = parse_file("../test/data/ravens/ravens_case3_lm_1230.json")
ravens_case3_lm_models = parse_file("../test/data/ravens/ravens_case3_lm_models.json")
ravens_case3_unbalanced_ZIPloads = parse_file("../test/data/ravens/ravens_case3_unbalanced_ZIPloads.json")

# 3-bus delta generators
ravens_case3_delta_gens = parse_file("../test/data/ravens/ravens_case3_delta_gens.json")

# 4-bus
ravens_case4_phase_drop = parse_file("../test/data/ravens/ravens_case4_phase_drop.json")

# 5-bus
ravens_case5_phase_drop = parse_file("../test/data/ravens/ravens_case5_phase_drop.json")

# 4-bus 2w transformer
ravens_ut_trans_2w_yy = parse_file("../test/data/ravens/ravens_ut_trans_2w_yy.json")
ravens_ut_trans_2w_dy_lead = parse_file("../test/data/ravens/ravens_ut_trans_2w_dy_lead.json")
ravens_ut_trans_2w_dy_lag = parse_file("../test/data/ravens/ravens_ut_trans_2w_dy_lag.json")
ravens_ut_trans_2w_dy_lead_small_series_impedance = parse_file("../test/data/ravens/ravens_ut_trans_2w_dy_lead_small_series_impedance.json")
ravens_ut_trans_2w_yy_bank = parse_file("../test/data/ravens/ravens_ut_trans_2w_yy_bank.json")
ravens_ut_trans_2w_yy_bank = parse_file("../test/data/ravens/ravens_ut_trans_2w_yy_bank.json")
ravens_ut_trans_2w_yy_unbanked = parse_file("../test/data/ravens/ravens_ut_trans_2w_yy_bank.json"; bank_transformers=false)
ravens_ut_trans_2w_yy_oltc = parse_file("../test/data/ravens/ravens_ut_trans_2w_yy_oltc.json")

# 4-bus 3w transformer
ravens_ut_trans_3w_dyy_1 = parse_file("../test/data/ravens/ravens_ut_trans_3w_dyy_1.json")
ravens_ut_trans_3w_dyy_2 = parse_file("../test/data/ravens/ravens_ut_trans_3w_dyy_2.json")
ravens_ut_trans_3w_dyy_3 = parse_file("../test/data/ravens/ravens_ut_trans_3w_dyy_3.json")
ravens_ut_trans_3w_dyy_3_loadloss = parse_file("../test/data/ravens/ravens_ut_trans_3w_dyy_3_loadloss.json")
ravens_trans_3w_center_tap = parse_file("../test/data/ravens/ravens_trans_3w_center_tap.json")
RLF_trans_3w_center_tap = parse_file("../test/data/ravens/RLF_trans_3w_center_tap.json")


# IEEE13
ravens_IEEE13_Assets = parse_file("../test/data/ravens/ravens_IEEE13_Assets.json")

ravens_IEEE13_RegControl = parse_file("../test/data/ravens/ravens_IEEE13_RegControl.json")
ravens_IEEE13_CapControl = parse_file("../test/data/ravens/ravens_IEEE13_CapControl.json")

# explicit neutral
ravens_test_gen_3ph_wye = parse_file("../test/data/ravens_en_validation_case_data/ravens_test_gen_3ph_wye.json")
ravens_test_switch = parse_file("../test/data/ravens_en_validation_case_data/ravens_test_switch.json")
ravens_test_gen_1ph_wye = parse_file("../test/data/ravens_en_validation_case_data/ravens_test_gen_1ph_wye.json")
ravens_test_trans_dy = parse_file("../test/data/ravens_en_validation_case_data/ravens_test_trans_dy.json")

# distribution transformer equivalent cases
ravens_dist_transformer = parse_file("../test/data/ravens/ravens_dist_transformer.json")

# IEEE13 Feeder
#DNC ravens_IEEE13_Feeder_engr = parse_file("../test/data/ravens/ravens_ieee13_feeder.json", multinetwork=true,time_series="daily", transformations=[remove_line_limits!, remove_transformer_limits!])

