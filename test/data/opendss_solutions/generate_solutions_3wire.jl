using Pkg
cd("test")
Pkg.activate("./")
# Pkg.add("OpenDSSDirect")
# Pkg.add("JSON")


"""
This script uses OpenDSSDirect to obtain voltage profiles for the validation test cases,
and saves them as json files which are used in the unit tests so that these do not have 
to add OpenDSSDirect as a dependency.
"""

import OpenDSSDirect; const ODD = OpenDSSDirect
using JSON

# relative path to data and solution folders
data_dir = "data/opendss"
solution_dir = "data/opendss_solutions"

"Uses OpenDSSDirect to obtain the voltage profile"
function get_soldss_opendssdirect(dss_path::AbstractString; tolerance=missing)
    dir = pwd()
    ODD.Basic.ClearAll()
    ODD.dss("compile $dss_path")

    if !ismissing(tolerance)
        ODD.Solution.Convergence(tolerance)
        ODD.Solution.Solve()
    end

    sol_dss = Dict{String, Any}()
    # buses
    sol_dss["bus"] = Dict{String, Any}()
    bnames = ODD.Circuit.AllBusNames()
    for bname in bnames
        ODD.Circuit.SetActiveBus(bname)
        sol_dss["bus"][bname] = Dict("vm"=>Dict{Int, Float64}(), "va"=>Dict{Int, Float64}())
        v = ODD.Bus.Voltages()
        for (i,c) in enumerate(ODD.Bus.Nodes())
            sol_dss["bus"][bname]["vm"][c] = abs(v[i])
            sol_dss["bus"][bname]["va"][c] = angle(v[i])
        end
    end
    # restore original directory
    cd(dir)
    return sol_dss
end


cases = ["case2_diag", "case3_balanced_basefreq", "case3_balanced_cap", "case3_balanced_isc", "case3_balanced_prop-order",
"case3_delta_gens", "case3_lm_models_2", "case3_unbalanced_1phase-pv", "case3_unbalanced_assym_swap", "case3_unbalanced_delta_loads", "case3_unbalanced_missingedge", 
"case3_unbalanced", "case4_phase_drop", "case5_phase_drop", "ut_trans_2w_dy_lag", "ut_trans_2w_dy_lead_small_series_impedance", "case3_unbalanced_switch",
"ut_trans_2w_dy_lead", "ut_trans_2w_yy_oltc", "ut_trans_2w_yy", "ut_trans_3w_dyy_1"]

cases_diff = setdiff(readdir(data_dir), cases.*".dss")

for case in cases
    sol_path = "$solution_dir/$case.json"
    # only create solutions that are not present yet
    # if !isfile(sol_path)
        data_path = "$data_dir/$case.dss"
        sol = get_soldss_opendssdirect(data_path, tolerance=1E-10)
        sol["dss_file"] = "$case.dss"
        open("$solution_dir/$case.json", "w") do f
            JSON.print(f, sol, 4)
        end
    # end
end


