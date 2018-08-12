var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#ThreePhasePowerModels.jl-Documentation-1",
    "page": "Home",
    "title": "ThreePhasePowerModels.jl Documentation",
    "category": "section",
    "text": "CurrentModule = ThreePhasePowerModels"
},

{
    "location": "index.html#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "ThreePowerModels.jl is a Julia/JuMP extension package to PowerModels.jl for modeling of Multi-Phase (with a focus on three-phase) power grids. "
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The latest stable release of PowerModels can be installed using the Julia package manager withPkg.add(\"ThreePhasePowerModels\")For the current development version, \"checkout\" this package withPkg.checkout(\"ThreePhasePowerModels\")At least one solver is required for running ThreePhasePowerModels.  The open-source solver Ipopt is recommended, as it is extremely fast, and can be used to solve a wide variety of the problems and network formulations provided in ThreePhasePowerModels.  The Ipopt solver can be installed via the package manager withPkg.add(\"Ipopt\")Test that the package works by runningPkg.test(\"ThreePhasePowerModels\")"
},

{
    "location": "quickguide.html#",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "page",
    "text": ""
},

{
    "location": "quickguide.html#Quick-Start-Guide-1",
    "page": "Getting Started",
    "title": "Quick Start Guide",
    "category": "section",
    "text": ""
},

{
    "location": "library.html#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "library.html#ThreePhasePowerModels.LPLinUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.LPLinUBFForm",
    "category": "type",
    "text": "LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current in scalar form\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPLinUBFPowerModel",
    "page": "Library",
    "title": "ThreePhasePowerModels.LPLinUBFPowerModel",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPLinUBFPowerModel-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.LPLinUBFPowerModel",
    "category": "method",
    "text": "default Lin3Distflow constructor for scalar form\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPdiagUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.LPdiagUBFForm",
    "category": "type",
    "text": "LinDist3Flow per Sankur et al 2016, using vector variables for power, voltage and current\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPdiagUBFPowerModel",
    "page": "Library",
    "title": "ThreePhasePowerModels.LPdiagUBFPowerModel",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPdiagUBFPowerModel-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.LPdiagUBFPowerModel",
    "category": "method",
    "text": "default LP unbalanced DistFlow constructor\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPfullUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.LPfullUBFForm",
    "category": "type",
    "text": "Simplified BFM per Gan and Low 2014, PSCC, using matrix variables for power, voltage and current\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPfullUBFPowerModel",
    "page": "Library",
    "title": "ThreePhasePowerModels.LPfullUBFPowerModel",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.LPfullUBFPowerModel-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.LPfullUBFPowerModel",
    "category": "method",
    "text": "default LP unbalanced DistFlow constructor\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SDPUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.SDPUBFForm",
    "category": "type",
    "text": "SDP BFM per Gan and Low 2014, PSCC\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SDPUBFPowerModel",
    "page": "Library",
    "title": "ThreePhasePowerModels.SDPUBFPowerModel",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SDPUBFPowerModel-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.SDPUBFPowerModel",
    "category": "method",
    "text": "default SDP unbalanced DistFlow constructor\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SOCConicUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.SOCConicUBFForm",
    "category": "type",
    "text": "SOC relaxation of SDPUBFForm per Kim, Kojima, & Yamashita 2003, cast as a SOC\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SOCConicUBFPowerModel",
    "page": "Library",
    "title": "ThreePhasePowerModels.SOCConicUBFPowerModel",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SOCConicUBFPowerModel-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.SOCConicUBFPowerModel",
    "category": "method",
    "text": "default SOC unbalanced DistFlow constructor\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SOCNLPUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.SOCNLPUBFForm",
    "category": "type",
    "text": "SOC relaxation of SDPUBFForm per Kim, Kojima, & Yamashita 2003, cast as an QCP\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SOCNLPUBFPowerModel",
    "page": "Library",
    "title": "ThreePhasePowerModels.SOCNLPUBFPowerModel",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.SOCNLPUBFPowerModel-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.SOCNLPUBFPowerModel",
    "category": "method",
    "text": "default SOC unbalanced DistFlow constructor\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_ac_tp_opf-Tuple{Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_ac_tp_opf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_ac_tp_pf-Tuple{Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_ac_tp_pf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_opf-Tuple{Dict{String,Any},Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_opf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_opf-Tuple{String,Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_opf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_opf_bf-Tuple{Dict{String,Any},Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_opf_bf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_opf_bf-Tuple{String,Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_opf_bf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_opf_pbs-Tuple{Dict{String,Any},Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_opf_pbs",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_opf_pbs-Tuple{String,Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_opf_pbs",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_ots-Tuple{Dict{String,Any},Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_ots",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_ots-Tuple{String,Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_ots",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_pf-Tuple{Dict{String,Any},Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_pf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_pf-Tuple{String,Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_pf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_pf_bf-Tuple{Dict{String,Any},Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_pf_bf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_pf_bf-Tuple{String,Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_pf_bf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_pf_pbs-Tuple{Dict{String,Any},Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_pf_pbs",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_tp_pf_pbs-Tuple{String,Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_tp_pf_pbs",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.AbstractConicUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.AbstractConicUBFForm",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.AbstractLPUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.AbstractLPUBFForm",
    "category": "type",
    "text": "Abstract form for linear unbalanced power flow models\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.AbstractNLPUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.AbstractNLPUBFForm",
    "category": "type",
    "text": "\n\n"
},

{
    "location": "library.html#PowerModels.constraint_branch_current-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPLinUBFForm",
    "page": "Library",
    "title": "PowerModels.constraint_branch_current",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_branch_current_setpoint-Tuple{Any,PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.add_branch_current_setpoint",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_branch_flow_setpoint-Tuple{Any,PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.add_branch_flow_setpoint",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_bus_slack_setpoint-Tuple{Any,PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.add_bus_slack_setpoint",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_bus_voltage_setpoint-Tuple{Any,PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.add_bus_voltage_setpoint",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_component!-Tuple{Dict,AbstractString,Dict}",
    "page": "Library",
    "title": "ThreePhasePowerModels.add_component!",
    "category": "method",
    "text": "add_component!(dss_data, ctype_name, compDict)\n\nAdds a component of type ctype_name with properties given by compDict to the existing dss_data structure. If a component of the same type has already been added to dss_data, the new component is appeneded to the existing array of components of that type, otherwise a new array is created.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_original_variables-Tuple{Any,PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.add_original_variables",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.add_property-Tuple{Dict,AbstractString,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.add_property",
    "category": "method",
    "text": "add_property(compDict, key, value)\n\nAdds a property to an existing component properties dictionary compDict given the key and value of the property. If a property of the same name already exists inside compDict, the original value is converted to an array, and the new value is appended to the end.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.adjust_sourcegen_bounds!-Tuple{Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.adjust_sourcegen_bounds!",
    "category": "method",
    "text": "adjust_sourcegen_bounds!(tppm_data)\n\nChanges the bounds for the sourcebus generator by checking the emergamps of all of the branches attached to the sourcebus and taking the sum of non-infinite values. Defaults to Inf if all emergamps connected to sourcebus are also Inf.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.assign_property!-Tuple{Dict,AbstractString,AbstractString,AbstractString,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.assign_property!",
    "category": "method",
    "text": "assign_property!(dss_data, cType, cName, propName, propValue)\n\nAssigns a property with name propName and value propValue to the component of type cType named cName in dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.calc_tp_voltage_product_bounds-Tuple{PowerModels.GenericPowerModel,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.calc_tp_voltage_product_bounds",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.check_duplicate_components!-Tuple{Dict}",
    "page": "Library",
    "title": "ThreePhasePowerModels.check_duplicate_components!",
    "category": "method",
    "text": "check_duplicate_components!(dss_data)\n\nFinds duplicate components in dss_data and merges up, meaning that older data (lower indices) is always overwritten by newer data (higher indices).\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.check_network_data-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.check_network_data",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_flow_losses-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPLinUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_flow_losses",
    "category": "method",
    "text": "Defines branch flow model power flow equations\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_kcl_shunt_slack-Tuple{PowerModels.GenericPowerModel,Int64}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_kcl_shunt_slack",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_kcl_shunt_slack-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:Union{PowerModels.AbstractACTForm, PowerModels.AbstractBFForm, PowerModels.AbstractWRForm, PowerModels.AbstractWRMForm}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_kcl_shunt_slack",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_kcl_shunt_slack-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Int64,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractACPForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_kcl_shunt_slack",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from-Tuple{PowerModels.GenericPowerModel,Int64}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractACPForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\np[f_idx] ==  (g+g_fr)/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))\nq[f_idx] == -(b+b_fr)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus]))\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractDCPForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\np[f_idx] == -b*(t[f_bus] - t[t_bus])\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:Union{PowerModels.AbstractACTForm, PowerModels.AbstractWRForm, PowerModels.AbstractWRMForm}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off-Tuple{PowerModels.GenericPowerModel,Int64}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractACPForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off",
    "category": "method",
    "text": "p[f_idx] == z*(g/tm*v[f_bus]^2 + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))\nq[f_idx] == z*(-(b+c/2)/tm*v[f_bus]^2 - (-b*tr-g*ti)/tm*(v[f_bus]*v[t_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr+b*ti)/tm*(v[f_bus]*v[t_bus]*sin(t[f_bus]-t[t_bus])))\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractDCPForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off",
    "category": "method",
    "text": "-b*(t[f_bus] - t[t_bus] + vad_min*(1-branch_z[i])) <= p[f_idx] <= -b*(t[f_bus] - t[t_bus] + vad_max*(1-branch_z[i]))\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractWRForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_from_on_off",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\np[f_idx] ==        g/tm*w_fr[i] + (-g*tr+b*ti)/tm*(wr[i]) + (-b*tr-g*ti)/tm*(wi[i])\nq[f_idx] == -(b+c/2)/tm*w_fr[i] - (-b*tr-g*ti)/tm*(wr[i]) + (-g*tr+b*ti)/tm*(wi[i])\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to-Tuple{PowerModels.GenericPowerModel,Int64}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractACPForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\np[t_idx] ==  (g+g_to)*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))\nq[t_idx] == -(b+b_to)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus]))\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractDCPForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to",
    "category": "method",
    "text": "Do nothing, this model is symmetric\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:Union{PowerModels.AbstractACTForm, PowerModels.AbstractWRForm, PowerModels.AbstractWRMForm}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off-Tuple{PowerModels.GenericPowerModel,Int64}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractACPForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off",
    "category": "method",
    "text": "p[t_idx] == z*(g*v[t_bus]^2 + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[t_bus]-t[f_bus])) + (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))\nq[t_idx] == z*(-(b+c/2)*v[t_bus]^2 - (-b*tr+g*ti)/tm*(v[t_bus]*v[f_bus]*cos(t[f_bus]-t[t_bus])) + (-g*tr-b*ti)/tm*(v[t_bus]*v[f_bus]*sin(t[t_bus]-t[f_bus])))\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractDCPForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off",
    "category": "method",
    "text": "Do nothing, this model is symmetric\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:PowerModels.AbstractWRForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_ohms_tp_yt_to_on_off",
    "category": "method",
    "text": "Creates Ohms constraints (yt post fix indicates that Y and T values are in rectangular form)\n\np[t_idx] ==        g*w_to[i] + (-g*tr-b*ti)/tm*(wr[i]) + (-b*tr+g*ti)/tm*(-wi[i])\nq[t_idx] == -(b+c/2)*w_to[i] - (-b*tr+g*ti)/tm*(wr[i]) + (-g*tr-b*ti)/tm*(-wi[i])\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_branch_current-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.AbstractLPUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_branch_current",
    "category": "method",
    "text": "Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_branch_current-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.SDPUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_branch_current",
    "category": "method",
    "text": "Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_branch_current-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any}, Tuple{T}} where T<:Union{ThreePhasePowerModels.SOCConicUBFForm, ThreePhasePowerModels.SOCNLPUBFForm}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_branch_current",
    "category": "method",
    "text": "Defines relationship between branch (series) power flow, branch (series) current and node voltage magnitude\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_flow_losses-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPLinUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_flow_losses",
    "category": "method",
    "text": "Defines branch flow model power flow loss equations\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_flow_losses-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPdiagUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_flow_losses",
    "category": "method",
    "text": "Defines branch flow model power flow equations\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_flow_losses-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPfullUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_flow_losses",
    "category": "method",
    "text": "Defines branch flow model power flow equations\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_flow_losses-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:Union{ThreePhasePowerModels.AbstractConicUBFForm, ThreePhasePowerModels.AbstractNLPUBFForm}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_flow_losses",
    "category": "method",
    "text": "Defines branch flow model power flow equations\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_theta_ref-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPdiagUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_theta_ref",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_theta_ref-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any}, Tuple{T}} where T<:Union{ThreePhasePowerModels.AbstractConicUBFForm, ThreePhasePowerModels.AbstractNLPUBFForm}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_theta_ref",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_theta_ref-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any}, Tuple{T}} where T<:Union{PowerModels.AbstractACPForm, PowerModels.AbstractACTForm, PowerModels.AbstractDCPForm}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_theta_ref",
    "category": "method",
    "text": "Creates phase angle constraints at reference buses\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_theta_ref-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any}, Tuple{T}} where T<:Union{PowerModels.AbstractACTForm, PowerModels.AbstractBFForm, PowerModels.AbstractWRForm, PowerModels.AbstractWRMForm}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_theta_ref",
    "category": "method",
    "text": "do nothing, no way to represent this in these variables\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_theta_ref-Union{Tuple{PowerModels.GenericPowerModel{T},Int64}, Tuple{T}} where T<:Union{ThreePhasePowerModels.AbstractConicUBFForm, ThreePhasePowerModels.AbstractNLPUBFForm}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_theta_ref",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_voltage",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64}, Tuple{T}} where T<:PowerModels.AbstractACPForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_voltage",
    "category": "method",
    "text": "do nothing, this model does not have complex voltage constraints\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64}, Tuple{T}} where T<:PowerModels.AbstractWRForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_voltage",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPLinUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference",
    "category": "method",
    "text": "Defines branch flow model power flow equations\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPdiagUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference",
    "category": "method",
    "text": "Defines voltage drop over a branch, linking from and to side voltage\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPfullUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference",
    "category": "method",
    "text": "Defines voltage drop over a branch, linking from and to side voltage\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:Union{ThreePhasePowerModels.AbstractConicUBFForm, ThreePhasePowerModels.AbstractNLPUBFForm}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_tp_voltage_magnitude_difference",
    "category": "method",
    "text": "Defines voltage drop over a branch, linking from and to side voltage\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_voltage_magnitude_difference-Tuple{PowerModels.GenericPowerModel,Int64}",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_voltage_magnitude_difference",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.constraint_voltage_magnitude_difference-Union{Tuple{PowerModels.GenericPowerModel{T},Int64,Int64,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.LPLinUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.constraint_voltage_magnitude_difference",
    "category": "method",
    "text": "Defines voltage drop over a branch, linking from and to side voltage magnitude\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createCapacitor",
    "page": "Library",
    "title": "ThreePhasePowerModels.createCapacitor",
    "category": "function",
    "text": "createCapacitor(bus1, name, bus2=0; kwargs)\n\nCreates a Dict{String,Any} containing all of the expected properties for a Capacitor. If bus2 is not specified, the capacitor will be treated as a shunt. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createGenerator-Tuple{Any,AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.createGenerator",
    "category": "method",
    "text": "createGenerator(bus1, name; kwargs...)\n\nCreates a Dict{String,Any} containing all of the expected properties for a Generator. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createLine-Tuple{Any,Any,AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.createLine",
    "category": "method",
    "text": "createLine(bus1, bus2, name; kwargs...)\n\nCreates a Dict{String,Any} containing all of the properties for a Line. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createLinecode-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.createLinecode",
    "category": "method",
    "text": "createLinecode(name; kwargs...)\n\nCreates a Dict{String,Any} containing all of the properties of a Linecode. See OpenDSS documentation for valid fields and ways to specify the different properties. DEPRECIATED: Calculation all done inside of createLine() due to Rg, Xg. Merge linecode values into line kwargs values BEFORE calling createLine(). This is now mainly used for parsing linecode dicts into correct data types.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createLoad-Tuple{Any,AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.createLoad",
    "category": "method",
    "text": "createLoad(bus1, name; kwargs...)\n\nCreates a Dict{String,Any} containing all of the expected properties for a Load. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createPVSystem-Tuple{Any,AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.createPVSystem",
    "category": "method",
    "text": "createPVSystem(bus1, name; kwargs...)\n\nCreates a Dict{String,Any} containing all of the expected properties for a PVSystem. See OpenDSS document https://github.com/tshort/OpenDSS/blob/master/Doc/OpenDSS%20PVSystem%20Model.doc for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createReactor",
    "page": "Library",
    "title": "ThreePhasePowerModels.createReactor",
    "category": "function",
    "text": "createReactor(bus1, name, bus2=0; kwargs...)\n\nCreates a Dict{String,Any} containing all of the expected properties for a Reactor. If bus2 is not specified Reactor is treated like a shunt. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createTransformer-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.createTransformer",
    "category": "method",
    "text": "createTransformer(name; kwargs...)\n\nCreates a Dict{String,Any} containing all of the expected properties for a Transformer. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.createVSource",
    "page": "Library",
    "title": "ThreePhasePowerModels.createVSource",
    "category": "function",
    "text": "createVSource(bus1, name, bus2=0; kwargs...)\n\nCreates a Dict{String,Any} containing all of the expected properties for a Voltage Source. If bus2 is not specified, VSource will be treated like a generator. Mostly used as sourcebus which represents the circuit. See OpenDSS documentation for valid fields and ways to specify the different properties.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.create_starbus-Tuple{Dict,Dict}",
    "page": "Library",
    "title": "ThreePhasePowerModels.create_starbus",
    "category": "method",
    "text": "creates a starbus from a 3-winding transformer\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.discover_buses-Tuple{Dict}",
    "page": "Library",
    "title": "ThreePhasePowerModels.discover_buses",
    "category": "method",
    "text": "discover_buses(dss_data)\n\nDiscovers all of the buses (not separately defined in OpenDSS), from \"lines\".\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_branch!-Tuple{Dict,Dict,Bool}",
    "page": "Library",
    "title": "ThreePhasePowerModels.dss2tppm_branch!",
    "category": "method",
    "text": "dss2tppm_branch!(tppm_data, dss_data)\n\nAdds PowerModels-style branches to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_bus!",
    "page": "Library",
    "title": "ThreePhasePowerModels.dss2tppm_bus!",
    "category": "function",
    "text": "dss2tppm_bus!(tppm_data, dss_data)\n\nAdds PowerModels-style buses to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_gen!-Tuple{Dict,Dict,Bool}",
    "page": "Library",
    "title": "ThreePhasePowerModels.dss2tppm_gen!",
    "category": "method",
    "text": "dss2tppm_gen!(tppm_data, dss_data)\n\nAdds PowerModels-style generators to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_load!-Tuple{Dict,Dict,Bool}",
    "page": "Library",
    "title": "ThreePhasePowerModels.dss2tppm_load!",
    "category": "method",
    "text": "dss2tppm_load!(tppm_data, dss_data)\n\nAdds PowerModels-style loads to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_pvsystem!-Tuple{Dict,Dict,Bool}",
    "page": "Library",
    "title": "ThreePhasePowerModels.dss2tppm_pvsystem!",
    "category": "method",
    "text": "dss2tppm_pvsystem!(tppm_data, dss_data)\n\nAdds PowerModels-style pvsystems to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_shunt!-Tuple{Dict,Dict,Bool}",
    "page": "Library",
    "title": "ThreePhasePowerModels.dss2tppm_shunt!",
    "category": "method",
    "text": "dss2tppm_shunt!(tppm_data, dss_data)\n\nAdds PowerModels-style shunts to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.dss2tppm_transformer!-Tuple{Dict,Dict,Bool}",
    "page": "Library",
    "title": "ThreePhasePowerModels.dss2tppm_transformer!",
    "category": "method",
    "text": "dss2tppm_transformer!(tppm_data, dss_data, import_all)\n\nAdds PowerModels-style transformers (branches) to tppm_data from dss_data.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.find_bus-Tuple{AbstractString,Dict}",
    "page": "Library",
    "title": "ThreePhasePowerModels.find_bus",
    "category": "method",
    "text": "find_bus(busname, tppm_data)\n\nFinds the index number of the bus in existing data from the given busname.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.find_component-Tuple{Dict,AbstractString,AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.find_component",
    "category": "method",
    "text": "find_component(tppm_data, name, compType)\n\nReturns the component of compType with name from data of type Dict{String,Array}.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.get_dtypes-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.get_dtypes",
    "category": "method",
    "text": "Returns a Dict{String,Type} for the desired component comp, giving all of the expected data types\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.get_linecode-Tuple{Dict,AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.get_linecode",
    "category": "method",
    "text": "returns the linecode with name id\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.get_pbs_solution-Tuple{PowerModels.GenericPowerModel,Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.get_pbs_solution",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.get_prop_name-Tuple{AbstractString,Int64}",
    "page": "Library",
    "title": "ThreePhasePowerModels.get_prop_name",
    "category": "method",
    "text": "get_prop_name(ctype, i)\n\nReturns the ith property name for a given component type ctype.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.get_prop_name-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.get_prop_name",
    "category": "method",
    "text": "get_prop_name(ctype)\n\nReturns the property names in order for a given component type ctype.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.get_solution_tp-Tuple{PowerModels.GenericPowerModel,Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.get_solution_tp",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.isa_array-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.isa_array",
    "category": "method",
    "text": "checks if data is an opendss-style array string\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.isa_conn-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.isa_conn",
    "category": "method",
    "text": "checks is a string is a connection by checking the values\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.isa_matrix-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.isa_matrix",
    "category": "method",
    "text": "checks if data is an opendss-style matrix string\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.isa_rpn-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.isa_rpn",
    "category": "method",
    "text": "detects if expr is Reverse Polish Notation expression\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.make_mpv!-Tuple{Dict{String,Any},String,Array{String,1}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.make_mpv!",
    "category": "method",
    "text": "collects several from_keys in an array and sets it to the to_key, removes from_keys\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.matlab_to_tppm-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.matlab_to_tppm",
    "category": "method",
    "text": "Converts a Matlab dict into a ThreePhasePowerModels dict\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.merge_dss!-Tuple{Dict{String,Array},Dict{String,Array}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.merge_dss!",
    "category": "method",
    "text": "merge_dss!(dss_prime, dss_to_add)\n\nMerges two (partially) parsed OpenDSS files to the same dictionary dss_prime. Used in cases where files are referenced via the \"compile\" or \"redirect\" OpenDSS commands inside the originating file.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.ml2pm_branch-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.ml2pm_branch",
    "category": "method",
    "text": "convert raw branch data into arrays\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.ml2pm_bus-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.ml2pm_bus",
    "category": "method",
    "text": "convert raw bus data into arrays\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.ml2pm_gen-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.ml2pm_gen",
    "category": "method",
    "text": "convert raw generator data into arrays\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.ml2pm_load-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.ml2pm_load",
    "category": "method",
    "text": "convert raw load data into arrays\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.ml2pm_shunt-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.ml2pm_shunt",
    "category": "method",
    "text": "convert raw shunt data into arrays\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.objective_min_bus_power_slack-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.objective_min_bus_power_slack",
    "category": "method",
    "text": "a quadratic penalty for bus power slack variables\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_array",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_array",
    "category": "function",
    "text": "parse matrices according to active nodes\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_array-Tuple{Type,AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_array",
    "category": "method",
    "text": "parse_array(dtype, data)\n\nParses a OpenDSS style array string data into a one dimensional array of type dtype. Array strings are capped by either brackets, single quotes, or double quotes, and elements are separated by spaces.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_buscoords-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_buscoords",
    "category": "method",
    "text": "parse_buscoords(file)\n\nParses a Bus Coordinate file, in either \"dat\" or \"csv\" formats, where in \"dat\", columns are separated by spaces, and in \"csv\" by commas. File expected to contain \"bus,x,y\" on each line.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_busname-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_busname",
    "category": "method",
    "text": "parse_busname(busname)\n\nParses busnames as defined in OpenDSS, e.g. \"primary.1.2.3.0\".\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_component",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_component",
    "category": "function",
    "text": "parse_component(component, properies, compDict=Dict{String,Any}())\n\nParses a component with properties into a compDict. If compDict is not defined, an empty dictionary will be used. Assumes that unnamed properties are given in order, but named properties can be given anywhere.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_conn-Tuple{String}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_conn",
    "category": "method",
    "text": "parses connection \"conn\" specification reducing to wye or delta\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_dss-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_dss",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_dss-Tuple{IOStream}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_dss",
    "category": "method",
    "text": "parse_dss(filename)\n\nParses a OpenDSS file given by filename into a Dict{Array{Dict}}. Only supports components and options, but not commands, e.g. \"plot\" or \"solve\". Will also parse files defined inside of the originating DSS file via the \"compile\", \"redirect\" or \"buscoords\" commands.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_dss_with_dtypes!",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_dss_with_dtypes!",
    "category": "function",
    "text": "parse_dss_with_dtypes!(dss_data, toParse)\n\nParses the data in keys defined by toParse in dss_data using types given by the default properties from the get_prop_default function.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_element_with_dtype-Tuple{Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_element_with_dtype",
    "category": "method",
    "text": "parses the raw dss values into their expected data types\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_file-Tuple{IOStream}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_file",
    "category": "method",
    "text": "parse_file(io)\n\nParses the IOStream of a file into a Three-Phase PowerModels data structure.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_file-Tuple{String}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_file",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_line",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_line",
    "category": "function",
    "text": "parse_line(elements, curCompDict=Dict{String,Any}())\n\nParses an already separated line given by elements (an array) of an OpenDSS file into curCompDict. If not defined, curCompDict is an empty dictionary.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_matlab-Tuple{IOStream}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_matlab",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_matlab-Tuple{String}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_matlab",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_matlab_string-Tuple{String}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_matlab_string",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_matrix-Tuple{Type,AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_matrix",
    "category": "method",
    "text": "parse_matrix(dtype, data)\n\nParses a OpenDSS style triangular matrix string data into a two dimensional array of type dtype. Matrix strings are capped by either parenthesis or brackets, rows are separated by \"|\", and columns are separated by spaces.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_matrix-Union{Tuple{Array{T,N} where N,Array{Bool,N} where N,Int64,Any}, Tuple{Array{T,N} where N,Array{Bool,N} where N,Int64}, Tuple{Array{T,N} where N,Array{Bool,N} where N}, Tuple{T}} where T",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_matrix",
    "category": "method",
    "text": "parse matrices according to active nodes\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_opendss-Tuple{Dict}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_opendss",
    "category": "method",
    "text": "Parses a Dict resulting from the parsing of a DSS file into a PowerModels usable format.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_opendss-Tuple{IOStream}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_opendss",
    "category": "method",
    "text": "Parses a DSS file into a PowerModels usable format.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_options-Tuple{Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_options",
    "category": "method",
    "text": "parse_options(options)\n\nParses options defined with the set command in OpenDSS.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_properties-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_properties",
    "category": "method",
    "text": "parse_properties(properties)\n\nParses a string of properties of a component type, character by character into an array with each element containing (if present) the property name, \"=\", and the property value.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.parse_rpn",
    "page": "Library",
    "title": "ThreePhasePowerModels.parse_rpn",
    "category": "function",
    "text": "parses Reverse Polish Notation expr\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_opf-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.post_tp_opf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_opf_bf-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.post_tp_opf_bf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_opf_pbs-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.post_tp_opf_pbs",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_ots-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.post_tp_ots",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_pf-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.post_tp_pf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_pf_bf-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.post_tp_pf_bf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.post_tp_pf_pbs-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.post_tp_pf_pbs",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.psd_to_psd_complex-Tuple{PowerModels.GenericPowerModel,Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.psd_to_psd_complex",
    "category": "method",
    "text": "complex SDP to SDP relaxation based on PSDness of principal minors\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.psd_to_psd_complex-Union{Tuple{PowerModels.GenericPowerModel{T},Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.AbstractConicUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.psd_to_psd_complex",
    "category": "method",
    "text": "complex SDP to SDP relaxation based on PSDness of principal minors, default is 3x3 SDP relaxation\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.psd_to_soc-Tuple{PowerModels.GenericPowerModel,Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.psd_to_soc",
    "category": "method",
    "text": "See section 4.3 in: Fazel, M., Hindi, H., & Boyd, S. P. (2001). A rank minimization heuristic with application to minimum order system approximation. Proc. American Control Conf., 6(2), 4734–4739. https://doi.org/10.1109/ACC.2001.945730\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.psd_to_soc-Tuple{PowerModels.GenericPowerModel,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.psd_to_soc",
    "category": "method",
    "text": "SDP to SOC relaxation of type 2 as described in: Kim, S., Kojima, M., & Yamashita, M. (2003). Second order cone programming relaxation of a positive semidefinite constraint. Optimization Methods and Software, 18(5), 535–541. https://doi.org/10.1080/1055678031000148696 Applied to real-value matrix\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.psd_to_soc-Union{Tuple{PowerModels.GenericPowerModel{T},Any}, Tuple{T}} where T<:ThreePhasePowerModels.AbstractConicUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.psd_to_soc",
    "category": "method",
    "text": "SDP to SOC relaxation of type 2 as described in: Kim, S., Kojima, M., & Yamashita, M. (2003). Second order cone programming relaxation of a positive semidefinite constraint. Optimization Methods and Software, 18(5), 535–541. https://doi.org/10.1080/1055678031000148696 Applied to real-value matrix\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.psd_to_soc_complex-Tuple{PowerModels.GenericPowerModel,Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.psd_to_soc_complex",
    "category": "method",
    "text": "SDP to SOC relaxation of type 2 as described in: Kim, S., Kojima, M., & Yamashita, M. (2003). Second order cone programming relaxation of a positive semidefinite constraint. Optimization Methods and Software, 18(5), 535–541. https://doi.org/10.1080/1055678031000148696\n\nApplied to complex-value matrix to obtain SOC: Andersen, M. S., Hansson, A., & Vandenberghe, L. (2014). Reduced-complexity semidefinite relaxations of optimal power flow problems. IEEE Trans. Power Syst., 29(4), 1855–1863.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.psd_to_soc_complex-Union{Tuple{PowerModels.GenericPowerModel{T},Any,Any}, Tuple{T}} where T<:ThreePhasePowerModels.AbstractConicUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.psd_to_soc_complex",
    "category": "method",
    "text": "SDP to SOC relaxation of type 2 as described in: Kim, S., Kojima, M., & Yamashita, M. (2003). Second order cone programming relaxation of a positive semidefinite constraint. Optimization Methods and Software, 18(5), 535–541. https://doi.org/10.1080/1055678031000148696\n\nApplied to complex-value matrix to obtain SOC: Andersen, M. S., Hansson, A., & Vandenberghe, L. (2014). Reduced-complexity semidefinite relaxations of optimal power flow problems. IEEE Trans. Power Syst., 29(4), 1855–1863.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.roll-Union{Tuple{Array{T,1},Int64}, Tuple{T}} where T<:Number",
    "page": "Library",
    "title": "ThreePhasePowerModels.roll",
    "category": "method",
    "text": "rolls a 1d array left or right by idx\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.run_dc_tp_pf-Tuple{Any,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.run_dc_tp_pf",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.set_default-Tuple{Dict{String,Any},String,Any}",
    "page": "Library",
    "title": "ThreePhasePowerModels.set_default",
    "category": "method",
    "text": "checks if the given dict has a value, if not, sets a default value\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.sqr-Tuple{Float64}",
    "page": "Library",
    "title": "ThreePhasePowerModels.sqr",
    "category": "method",
    "text": "Squares x, for parsing Reverse Polish Notation\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.strip_comments-Tuple{AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.strip_comments",
    "category": "method",
    "text": "Strips comments, defined by \"!\" from the ends of lines\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.strip_lines-Tuple{Array}",
    "page": "Library",
    "title": "ThreePhasePowerModels.strip_lines",
    "category": "method",
    "text": "strips lines that are either commented (block or single) or empty\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.to_sym_keys-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.to_sym_keys",
    "category": "method",
    "text": "converts Dict{String,Any} to Dict{Symbol,Any} for passing as kwargs\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.translate_version!-Tuple{Dict{String,Any}}",
    "page": "Library",
    "title": "ThreePhasePowerModels.translate_version!",
    "category": "method",
    "text": "Translates legacy versions into current version format\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_active_bus_power_slack-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_active_bus_power_slack",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_bus_power_slack-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_bus_power_slack",
    "category": "method",
    "text": "generates variables for both active and reactive slack at each bus\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_lower_triangle_active_branch_flow-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_lower_triangle_active_branch_flow",
    "category": "method",
    "text": "variable: p_lt[l,i,j] for (l,i,j) in arcs\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_lower_triangle_reactive_branch_flow-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_lower_triangle_reactive_branch_flow",
    "category": "method",
    "text": "variable: q_lt[l,i,j] for (l,i,j) in arcs\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_reactive_bus_power_slack-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_reactive_bus_power_slack",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_branch_flow-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_tp_branch_flow",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_branch_flow-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:ThreePhasePowerModels.LPdiagUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_tp_branch_flow",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_branch_flow-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:ThreePhasePowerModels.LPfullUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_tp_branch_flow",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_branch_flow-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:Union{ThreePhasePowerModels.AbstractConicUBFForm, ThreePhasePowerModels.AbstractNLPUBFForm}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_tp_branch_flow",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_tp_voltage",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:PowerModels.AbstractWRForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_tp_voltage",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:ThreePhasePowerModels.LPLinUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_tp_voltage",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage_magnitude_sqr-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_tp_voltage_magnitude_sqr",
    "category": "method",
    "text": "variable: w[i] >= 0 for i in buses\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage_prod_hermitian-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:ThreePhasePowerModels.LPdiagUBFForm",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_tp_voltage_prod_hermitian",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage_prod_hermitian-Union{Tuple{PowerModels.GenericPowerModel{T}}, Tuple{T}} where T<:Union{ThreePhasePowerModels.AbstractConicUBFForm, ThreePhasePowerModels.AbstractNLPUBFForm}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_tp_voltage_prod_hermitian",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_tp_voltage_product-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_tp_voltage_product",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_upper_triangle_active_branch_flow-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_upper_triangle_active_branch_flow",
    "category": "method",
    "text": "variable: p_ut[l,i,j] for (l,i,j) in arcs\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.variable_upper_triangle_reactive_branch_flow-Tuple{PowerModels.GenericPowerModel}",
    "page": "Library",
    "title": "ThreePhasePowerModels.variable_upper_triangle_reactive_branch_flow",
    "category": "method",
    "text": "variable: q_ut[l,i,j] for (l,i,j) in arcs\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.where_is_comp-Tuple{Array,AbstractString}",
    "page": "Library",
    "title": "ThreePhasePowerModels.where_is_comp",
    "category": "method",
    "text": "where_is_comp(data, comp_id)\n\nFinds existing component of id comp_id in array of data and returns index. Assumes all components in data are unique.\n\n\n\n"
},

{
    "location": "library.html#ThreePhasePowerModels.jl-Library-1",
    "page": "Library",
    "title": "ThreePhasePowerModels.jl Library",
    "category": "section",
    "text": "Modules = [ThreePhasePowerModels]"
},

{
    "location": "developer.html#",
    "page": "Developer",
    "title": "Developer",
    "category": "page",
    "text": ""
},

{
    "location": "developer.html#Developer-Documentation-1",
    "page": "Developer",
    "title": "Developer Documentation",
    "category": "section",
    "text": ""
},

]}
