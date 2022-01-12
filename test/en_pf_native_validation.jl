@info "running explicit neutral power flow tests with native julia power flow solver"

function vsource_correction!(data_eng)
    data_eng["voltage_source"]["source"]["rs"][4,4] = data_eng["voltage_source"]["source"]["rs"][1,1]
    data_eng["voltage_source"]["source"]["rs"][1:3,4] .= data_eng["voltage_source"]["source"]["rs"][1,2]
    data_eng["voltage_source"]["source"]["rs"][4,1:3] .= data_eng["voltage_source"]["source"]["rs"][1,2]
    data_eng["voltage_source"]["source"]["xs"][4,4] = data_eng["voltage_source"]["source"]["xs"][1,1]
    data_eng["voltage_source"]["source"]["xs"][1:3,4] .= data_eng["voltage_source"]["source"]["xs"][1,2]
    data_eng["voltage_source"]["source"]["xs"][4,1:3] .= data_eng["voltage_source"]["source"]["xs"][1,2]
    return nothing
end

# point to data and solution directory
data_dir = "data/en_validation_case_data"
solution_dir = "data/en_validation_case_solutions"
# infer cases from files defined in data dir
cases = [x[1:end-4] for x in readdir(data_dir) if endswith(x, ".dss")]



@testset "en pf native opendss validation" begin

    for (case_idx,case) in enumerate(cases)

        @testset "case $case" begin
            case_path = "$data_dir/$case.dss"

            # transformations = haskey(case_transformations, case) ? case_transformations[case] : [remove_all_bounds!]
            data_eng = parse_file(case_path, transformations=[transform_loops!])
            vsource_correction!(data_eng)

            data_math = transform_data_model(data_eng;kron_reduce=false)

            add_start_vrvi!(data_math)
            v_start = PMD._bts_to_start_voltage(data_math)

            res = compute_pf(data_math)

            # obtain solution from dss
            sol_dss = open("$solution_dir/$case.json", "r") do f
                JSON.parse(f)
            end

            sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
            # @assert res["termination_status"]==LOCALLY_SOLVED

            v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false)
            @test v_maxerr_pu <= 1E-6

        end
    end
end

