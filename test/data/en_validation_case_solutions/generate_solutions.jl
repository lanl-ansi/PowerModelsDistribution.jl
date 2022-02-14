"""
This script uses OpenDSSDirect to obtain voltage profiles for the validation test cases,
and saves them as json files which are used in the unit tests so that these do not have 
to add OpenDSSDirect as a dependency.

The script assumes pwd() = test/data

If the test case is a multiperiod one, the dir has to be set to the load_profile.csv dir.
"""

using Pkg
cd("test/data")
Pkg.activate("./")
Pkg.add("OpenDSSDirect")
Pkg.add("JSON")

import OpenDSSDirect; const ODD = OpenDSSDirect
using JSON

# relative path to data and solution folders
data_dir = "en_validation_case_data"
solution_dir = "en_validation_case_solutions"

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

cases = [x[1:end-4] for x in readdir(data_dir) if endswith(x, ".dss")]

for case in cases
    sol_path = "$solution_dir/$case.json"
    # only create solutions that are not present yet
    if !isfile(sol_path)
        data_path = "$data_dir/$case.dss"
        sol = get_soldss_opendssdirect(data_path, tolerance=1E-10)
        sol["dss_file"] = "$case.dss"
        open("$solution_dir/$case.json", "w") do f
            JSON.print(f, sol, 4)
        end
    end
end

