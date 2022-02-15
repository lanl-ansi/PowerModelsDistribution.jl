@info "running explicit neutral power flow tests with native julia power flow solver"
using Pkg
cd("test")
Pkg.activate("./")
# Pkg.add("../#four-wire-native-pf")
# Pkg.add("Test")
# Pkg.add("JSON")
# Pkg.add("Ipopt")
using PowerModelsDistribution
using JSON
using Ipopt
using Test

function vsource_correction!(data_eng)
    if haskey(data_eng, "multinetwork")# && data_eng["multinetwork"]
        for (n,nw) in data_eng["nw"]
            nw["voltage_source"]["source"]["rs"][4,4] = nw["voltage_source"]["source"]["rs"][1,1]
            nw["voltage_source"]["source"]["rs"][1:3,4] .= nw["voltage_source"]["source"]["rs"][1,2]
            nw["voltage_source"]["source"]["rs"][4,1:3] .= nw["voltage_source"]["source"]["rs"][1,2]
            nw["voltage_source"]["source"]["xs"][4,4] = nw["voltage_source"]["source"]["xs"][1,1]
            nw["voltage_source"]["source"]["xs"][1:3,4] .= nw["voltage_source"]["source"]["xs"][1,2]
            nw["voltage_source"]["source"]["xs"][4,1:3] .= nw["voltage_source"]["source"]["xs"][1,2]
        end
    else
        data_eng["voltage_source"]["source"]["rs"][4,4] = data_eng["voltage_source"]["source"]["rs"][1,1]
        data_eng["voltage_source"]["source"]["rs"][1:3,4] .= data_eng["voltage_source"]["source"]["rs"][1,2]
        data_eng["voltage_source"]["source"]["rs"][4,1:3] .= data_eng["voltage_source"]["source"]["rs"][1,2]
        data_eng["voltage_source"]["source"]["xs"][4,4] = data_eng["voltage_source"]["source"]["xs"][1,1]
        data_eng["voltage_source"]["source"]["xs"][1:3,4] .= data_eng["voltage_source"]["source"]["xs"][1,2]
        data_eng["voltage_source"]["source"]["xs"][4,1:3] .= data_eng["voltage_source"]["source"]["xs"][1,2]        
    end
    return nothing
end


function multinetwork_data_math_correction!(data_math::Dict{String, Any})
    @assert data_math["multinetwork"]
    @assert data_math["data_model"]==MATHEMATICAL
    for (nw, dm) in data_math["nw"]
        dm["data_model"] = MATHEMATICAL
        dm["map"] = data_math["map"]
        dm["bus_lookup"] = data_math["bus_lookup"]
    end
    return nothing
end


# point to data and solution directory
data_dir = "data/en_validation_case_data"
solution_dir = "data/en_validation_case_solutions"
# infer cases from files defined in data dir
cases = [x[1:end-4] for x in readdir(data_dir) if endswith(x, ".dss")]
filter!(e->e≠"test_line_6w", cases)

@testset "en pf native opendss validation four wire" begin

    for (case_idx,case) in enumerate(cases)

        @testset "case $case" begin
            case_path = "$data_dir/$case.dss"

            data_eng = parse_file(case_path, transformations=[transform_loops!])
            vsource_correction!(data_eng)

            data_math = transform_data_model(data_eng;kron_reduce=false)

            res = compute_pf(data_math)

            # obtain solution from dss
            sol_dss = open("$solution_dir/$case.json", "r") do f
                JSON.parse(f)
            end

            sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
            @test res["termination_status"]==CONVERGED

            v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false)
            @test v_maxerr_pu <= 1E-8

        end
    end
end


@testset "en pf native opendss validation for max iteration limit" begin
    case_path = "$data_dir/test_load_3ph_wye_cp.dss"

    data_eng = parse_file(case_path, transformations=[transform_loops!])
    vsource_correction!(data_eng)

    data_math = transform_data_model(data_eng;kron_reduce=false)

    res = compute_pf(data_math; max_iter=5)
    @test res["termination_status"] == ITERATION_LIMIT
end


@testset "en pf native opendss validation for 6 wire test case" begin
    case_path = "$data_dir/test_line_6w.dss"

    data_eng = parse_file(case_path, transformations=[transform_loops!])
    vsource_correction!(data_eng)

    data_math = transform_data_model(data_eng;kron_reduce=false)

    try
        res = compute_pf(data_math)
        @test false
    catch Error
        @test true
    end

end


@testset "en pf native opendss validation for multinetwork test case" begin
    case_path = "$data_dir/opendss/case3_balanced.dss"
    case = "case3_balanced"
    # ipopt_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
    # case3_balanced = parse_file("data/opendss/case3_balanced.dss")

    eng_ts = make_multinetwork(case3_balanced)
    vsource_correction!(eng_ts)

    ## This section is to validate that the new feature does not break multinetwork 
    result_mn = solve_mn_mc_opf(eng_ts, ACPUPowerModel, ipopt_solver)

    @test result_mn["termination_status"] == LOCALLY_SOLVED
    @test isapprox(result_mn["solution"]["nw"]["1"]["bus"]["loadbus"]["vm"][1], 0.229971669, atol=1E-6)
    @test isapprox(result_mn["solution"]["nw"]["10"]["bus"]["loadbus"]["vm"][1], 0.22977482, atol=1E-6)


    ## This section is to validate that the native pf supports multinetwork
    data_math = transform_data_model(eng_ts;kron_reduce=false)
    multinetwork_data_math_correction!(data_math)   # this is done to avoide error in pf.jl line 378,
                                                    # however, it seems that compute_pf is not coded for multiperiod from the beginning?

    res = compute_pf(data_math)
    @test res["termination_status"]["10"]==CONVERGED

    sol_dss = open("$solution_dir/$case.json", "r") do f
        JSON.parse(f)
    end
    
    sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
    
    # data_math["nw"]["10"]["bus_lookup"] = data_math["bus_lookup"]["10"]
    v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd["nw"]["10"], eng_ts["nw"]["10"], data_math["nw"]["10"], verbose=false)
    @test v_maxerr_pu <= 1E-8

end


##
# point to data and solution directory
data_dir = "data/opendss"
solution_dir = "data/opendss_solutions"
# infer cases from files defined in data dir
cases = [x[1:end-5] for x in readdir(solution_dir) if endswith(x, ".json")]

@testset "en pf native opendss validation three wire" begin

    for (case_idx,case) in enumerate(cases)

        @testset "case $case" begin
            case = "case3_lm_models"
            case = "case3_delta_gens"
            case_path = "$data_dir/$case.dss"

            data_eng = parse_file(case_path, transformations=[transform_loops!])
            vsource_correction!(data_eng)

            data_math = transform_data_model(data_eng;kron_reduce=false)

            res = compute_pf(data_math)

            # obtain solution from dss
            sol_dss = open("$solution_dir/$case.json", "r") do f
                JSON.parse(f)
            end

            sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
            @test res["termination_status"]==CONVERGED

            v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false)
            @test v_maxerr_pu <= 1E-8

        end
    end
end
