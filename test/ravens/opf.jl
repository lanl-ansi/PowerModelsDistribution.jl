@info "running optimal power flow (opf) tests"

@testset "test opf" begin
    @testset "4-bus phase drop acp opf" begin
        math = transform_data_model(ravens_case4_phase_drop)
        result = solve_mc_opf(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)
        
        @test all(isapprox(sum(result["voltage_source"]["source"]["pg"]), 18.2632; atol=10*1e-2))
        @test all(isapprox(sum(result["voltage_source"]["source"]["qg"]), 9.02334; atol=10*1e-2))

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["loadbus1"]["vm"] ./ vbase, 0.98995; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["loadbus1"]["va"], 0.27; atol=10*1e-2))

        @test all(isapprox.(result["bus"]["loadbus2"]["vm"] ./ vbase, 0.98803; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["loadbus2"]["va"], -119.74; atol=10*1e-2))

        @test all(isapprox.(result["bus"]["loadbus3"]["vm"] ./ vbase, 0.98611; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["loadbus3"]["va"], 120.25; atol=10*1e-2))
    end

    @testset "4-bus phase drop acr opf" begin
        math = transform_data_model(ravens_case4_phase_drop)
        result = solve_mc_opf(math, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test all(isapprox(sum(result["voltage_source"]["source"]["pg"]), 18.2632; atol=10*1e-2))
        @test all(isapprox(sum(result["voltage_source"]["source"]["qg"]), 9.02334; atol=10*1e-2))

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["loadbus1"]["vm"] ./ vbase, 0.98995; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["loadbus1"]["va"], 0.27; atol=10*1e-2))

        @test all(isapprox.(result["bus"]["loadbus2"]["vm"] ./ vbase, 0.98803; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["loadbus2"]["va"], -119.74; atol=10*1e-2))

        @test all(isapprox.(result["bus"]["loadbus3"]["vm"] ./ vbase, 0.98611; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["loadbus3"]["va"], 120.25; atol=10*1e-2))
    end

    @testset "5-bus phase drop acp opf" begin
        math = transform_data_model(ravens_case5_phase_drop)
        result = solve_mc_opf(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test all(isapprox(sum(result["voltage_source"]["source"]["pg"]),  59.9363; atol=10*1e-2))
        @test all(isapprox(sum(result["voltage_source"]["source"]["qg"]), -33.5395; atol=10*1e-2))

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["midbus"]["vm"] ./ vbase, [0.97351, 0.96490, 0.95646]; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["midbus"]["va"], [-1.26, -121.31, 118.17]; atol=10*1e-2))
    end

    @testset "5-bus phase drop acr opf" begin
        math = transform_data_model(ravens_case5_phase_drop)
        result = solve_mc_opf(math, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test all(isapprox(sum(result["voltage_source"]["source"]["pg"]),  59.9363; atol=10*1e-2))
        @test all(isapprox(sum(result["voltage_source"]["source"]["qg"]), -33.5395; atol=10*1e-2))

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["midbus"]["vm"] ./ vbase, [0.97351, 0.96490, 0.95646]; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["midbus"]["va"], [-1.26, -121.31, 118.17]; atol=10*1e-2))
    end

    @testset "5-bus phase drop dcp opf" begin
        math = transform_data_model(ravens_case5_phase_drop)
        result = solve_mc_opf(math, DCPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test all(isapprox(sum(result["voltage_source"]["source"]["pg"]), 54.0; atol=10*1e-1))
        @test all(isapprox.(result["bus"]["midbus"]["va"], [-15.19, -136.64, 93.24]; atol=10*1e-1))
    end

    @testset "5-bus phase drop nfa opf" begin
        math = transform_data_model(ravens_case5_phase_drop)
        result = solve_mc_opf(math, NFAUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test all(isapprox(sum(result["voltage_source"]["source"]["pg"]), 54.0; atol=10*1e-1))
    end

    # @testset "test warm start from dss voltages export" begin
    #     eng = deepcopy(case5_phase_drop)
    #     math = transform_data_model(ravens_case5_phase_drop)
    #     add_voltage_starts!(eng, "../test/data/opendss/case5_voltage.csv")
    #     apply_voltage_bounds_math!(eng)

    #     for (_,bus) in eng["bus"]
    #         nt = length(bus["terminals"])
    #         @test haskey(bus, "vm_start") && length(bus["vm_start"]) == nt
    #         @test haskey(bus, "va_start") && length(bus["va_start"]) == nt
    #     end

    #     math = transform_data_model(eng)

    #     for (_,bus) in math["bus"]
    #         if startswith(bus["source_id"], "bus")
    #             @test bus["vm_start"] <= bus["vmax"] && bus["vm_start"] >= bus["vmin"]
    #         end
    #     end

    #     pm = instantiate_mc_model(math, ACPUPowerModel, build_mc_opf)
    #     for (i,vm) in var(pm, nw_id_default, :vm)
    #         if startswith(math["bus"]["$i"]["source_id"], "bus")
    #             @test all(JuMP.start_value(vm[t]) == math["bus"]["$i"]["vm_start"][idx] for (idx,t) in enumerate(math["bus"]["$i"]["terminals"]))
    #         end
    #     end

    #     for (i,va) in var(pm, nw_id_default, :va)
    #         if startswith(math["bus"]["$i"]["source_id"], "bus")
    #             @test all(JuMP.start_value(va[t]) == math["bus"]["$i"]["va_start"][idx] for (idx,t) in enumerate(math["bus"]["$i"]["terminals"]))
    #         end
    #     end

    #     result = solve_mc_opf(eng, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

    #     vbase = case5_phase_drop["settings"]["vbases_default"]["sourcebus"]
    #     @test all(isapprox.(result["bus"]["midbus"]["vm"] ./ vbase, [0.97352, 0.9649, 0.95646]; atol=10*1e-4))
    #     @test all(isapprox.(result["bus"]["midbus"]["va"], [-1.3, -121.3, 118.2]; atol=10*1e-1))
    # end

    @testset "2-bus diagonal acp opf" begin
        math = transform_data_model(ravens_case2_diag)
        result = solve_mc_opf(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 18.209; atol=10*1e-2)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]),  0.208; atol=10*1e-2)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["primary"]["vm"] ./ vbase, 0.984406; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["primary"]["va"], [-0.66, -120.66, 119.34]; atol=10*1e-2))
    end

    @testset "3-bus balanced acp opf" begin
        math = transform_data_model(ravens_case3_balanced)
        result = solve_mc_opf(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 18.345; atol=10*1e-2)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]),  9.194; atol=10*1e-1)

        @test all(isapprox.(result["bus"]["sourcebus"]["vm"], 0.229993; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["sourcebus"]["va"], [0.0, -120.0, 120.0]; atol=10*1e-2))

        @test all(isapprox.(result["bus"]["primary"]["vm"], 0.227932; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["primary"]["va"], [-0.08, -120.08, 119.92]; atol=10*0.2))

        @test all(isapprox.(result["bus"]["loadbus"]["vm"], 0.225537; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["loadbus"]["va"], [-0.17, -120.17, 119.83]; atol=10*0.2))
    end

    @testset "3-bus unbalanced acp opf" begin
        math = transform_data_model(ravens_case3_unbalanced)
        result = solve_mc_opf(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["sourcebus"]["vm"] ./ vbase, [0.9959, 0.9959, 0.9959]; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["sourcebus"]["va"], [0.0, -120.0, 120.0]; atol=10*1e-2))

        @test all(isapprox.(result["bus"]["primary"]["vm"] ./ vbase, [0.98094, 0.989365, 0.987043]; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["primary"]["va"], [-0.22, -120.11, 120.12]; atol=10*1e-2))

        @test all(isapprox.(result["bus"]["loadbus"]["vm"] ./ vbase, [0.96355, 0.981767, 0.976786]; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["loadbus"]["va"], [-0.48, -120.24, 120.27]; atol=10*1e-2))

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 21.4812; atol=10*1e-2)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), 9.27263; atol=10*1e-2)
    end

    @testset "3-bus balanced isc acp opf" begin
        math = transform_data_model(ravens_case3_balanced_isc)
        result = solve_mc_opf(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 18.345; atol=10*1e-2)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]),  9.194; atol=10*1e-2)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["primary"]["vm"] ./ vbase, 0.986953; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["primary"]["va"], [-0.03, -120.03, 119.97]; atol=10*1e-2))

        @test all(isapprox.(result["bus"]["loadbus"]["vm"] ./ vbase, 0.976586; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["loadbus"]["va"], [-0.07, -120.07, 119.93]; atol=10*1e-2))
    end

    @testset "3-bus balanced pv acp opf" begin
        math = transform_data_model(ravens_case3_balanced_pv)
        result = solve_mc_opf(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test sum(result["voltage_source"]["source"]["pg"]) < 0.0
        @test sum(result["voltage_source"]["source"]["qg"]) < 5.0

        @test isapprox(sum(result["solar"]["pv1"]["pg"]), 18.374; atol=10*1e-1)
        @test isapprox(sum(result["solar"]["pv1"]["qg"]),  4.994; atol=10*1e-0)

        # TODO improve accuracy of pv model
        # vbase = case3_balanced_pv["settings"]["vbases_default"]["sourcebus"]
        # @test all(isapprox.(result["bus"]["primary"]["vm"] ./ vbase, 0.996456; atol=10*1e-4))
    end

    @testset "3-bus unbalanced single-phase pv acp opf" begin
        math = transform_data_model(ravens_case3_unbalanced_1phase_pv)
        result = solve_mc_opf(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 18.38728; atol=10*1e-1)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]),  9.21713; atol=10*1e-0)

        @test all(isapprox(sum(result["solar"]["pv1"]["pg"]), 1.9947; atol=10*1e-1))
        @test all(isapprox(sum(result["solar"]["pv1"]["qg"]), 0.0; atol=10*2e-0))

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox(result["bus"]["primary"]["vm"] ./ vbase, [0.990984, 0.991149, 0.991134]; atol=10*1e-3))
        @test all(isapprox(result["bus"]["primary"]["va"], [-0.03, -120.03, 119.97]; atol=10*2e-1))
    end

    @testset "3-bus balanced capacitor acp opf" begin
        math = transform_data_model(ravens_case3_balanced_cap)
        result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        
        @test all(isapprox.(result["bus"]["primary"]["vm"] ./ vbase, 0.99127; atol=10*1e-4))
        @test all(isapprox.(result["bus"]["loadbus"]["vm"] ./ vbase, 0.98588; atol=10*1e-4))
    end

    @testset "3w dyy transformer nfa opf" begin
        math = transform_data_model(ravens_ut_trans_3w_dyy_1)
        result = solve_mc_opf(math, NFAUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 616.0; atol=10*1.0)
    end

    @testset "vm and va start with ivr and acr opf" begin
        pmd_eng = deepcopy(case3_unbalanced)
        pmd_math = transform_data_model(pmd_eng)

        sol_ivr = solve_mc_opf(pmd_math, IVRUPowerModel, ipopt_solver)
        sol_acr = solve_mc_opf(pmd_math, ACRUPowerModel, ipopt_solver)

        pmd_math["bus"]["4"]["vm_start"] = pmd_math["bus"]["4"]["vm"]
        pmd_math["bus"]["4"]["va_start"] = pmd_math["bus"]["4"]["va"]
        pmd_math["bus"]["2"]["vm_start"] = [0.9959, 0.9959, 0.9959]
        pmd_math["bus"]["2"]["va_start"] = [0.00, -2.0944, 2.0944]

        sol_ivr_with_start = solve_mc_opf(pmd_math, IVRUPowerModel, ipopt_solver)
        sol_acr_with_start = solve_mc_opf(pmd_math, ACRUPowerModel, ipopt_solver)

        @test isapprox(sol_ivr["objective"], sol_ivr_with_start["objective"]; atol=10*1e-5)
        @test isapprox(sol_acr["objective"], sol_acr_with_start["objective"]; atol=10*1e-5)
    end

    @testset "assign start value per connection iteration" begin
        data = transform_data_model(case3_unbalanced)

        #add a single-phase generator and assign a single-phase start value to it
        data["gen"]["2"] = Dict{String, Any}(
                            "pg_start"      => [-0.012],
                            "qg_start"      => [-0.006],
                            "model"         => 2,
                            "connections"   => [2],
                            "shutdown"      => 0.0,
                            "startup"       => 0.0,
                            "configuration" => WYE,
                            "name"          => "single_ph_generator",
                            "gen_bus"       => 3,
                            "pmax"          => [Inf],
                            "vbase"         => 0.23094,
                            "index"         => 2,
                            "cost"          => [0.5, 0.0],
                            "gen_status"    => 1,
                            "qmax"          => [Inf],
                            "qmin"          => [-Inf],
                            "pmin"          => [-Inf],
                            "ncost"         => 2
                            )
        result = solve_mc_opf(data, ACPUPowerModel, optimizer_with_attributes(Ipopt.Optimizer, "sb"=>"yes","print_level"=>0))

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
    end

    @testset "3-bus unbalanced test w/ gen nl costs above quadratic" begin
        data = transform_data_model(case3_unbalanced)

        #add a single-phase generator and assign a single-phase start value to it
        data["gen"]["2"] = Dict{String, Any}(
                            "pg_start"      => [-0.012],
                            "qg_start"      => [-0.006],
                            "model"         => 2,
                            "connections"   => [2],
                            "shutdown"      => 0.0,
                            "startup"       => 0.0,
                            "configuration" => WYE,
                            "name"          => "single_ph_generator",
                            "gen_bus"       => 3,
                            "pmax"          => [Inf],
                            "vbase"         => 0.23094,
                            "index"         => 2,
                            "cost"          => [0.2, 0.5, 1.1, 1.4, 1.0],
                            "gen_status"    => 1,
                            "qmax"          => [Inf],
                            "qmin"          => [-Inf],
                            "pmin"          => [-Inf],
                            "ncost"         => 5
                            )
        result = solve_mc_opf(data, ACPUPowerModel, ipopt_solver)
        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.9666; atol = 1e-4)
    end

    @testset "3-bus unbalanced fotp opf with yy transformer" begin
        math = transform_data_model(ravens_ut_trans_2w_yy)
        result = solve_mc_opf(math, FOTPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == LOCALLY_INFEASIBLE || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 467.547; atol=10*40)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), 484.327; atol=10*30)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["1"]["vm"] ./ vbase, [0.87451, 0.8613, 0.85348]; atol=10*2e-2))
        @test all(isapprox.(result["bus"]["1"]["va"], [-0.1, -120.4, 119.8]; atol=10*2e-1))
    end

    @testset "3-bus unbalanced fotp opf with dy transformer" begin
        math = transform_data_model(ravens_ut_trans_2w_dy_lag)
        result = solve_mc_opf(math, FOTPUPowerModel, ipopt_solver)
        
        @test result["termination_status"] == LOCALLY_INFEASIBLE || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 467.699; atol=10*200)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), 485.553; atol=10*60)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]


        @test all(isapprox.(result["bus"]["1"]["vm"] ./ vbase, [0.92092, 0.91012, 0.90059]; atol=10*3e-1))
    end

    @testset "3-bus unbalanced fotp opf with voltage-dependent loads" begin
        math = transform_data_model(ravens_case3_unbalanced_delta_loads)
        result = solve_mc_opf(math, FOTPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_INFEASIBLE || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"])*2, 42.0464; atol=10*1)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), 18.1928; atol=10*1)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["loadbus"]["vm"] ./ vbase, [0.94105, 0.95942, 0.95876]; atol=10*1e-2))
        @test all(isapprox.(result["bus"]["loadbus"]["va"], [-0.9, -120.3, 120.2]; atol=10*5e-1))
    end

    @testset "3-bus unbalanced fotr opf with yy transformer" begin
        math = transform_data_model(ravens_ut_trans_2w_yy)
        result = solve_mc_opf(math, FOTRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 467.547; atol=10*70)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), 484.327; atol=10*90)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["1"]["vm"] ./ vbase, [0.87451, 0.8613, 0.85348]; atol=10*2e-2))
        @test all(isapprox.(result["bus"]["1"]["va"], [-0.1, -120.4, 119.8]; atol=10*2e-1))
    end

    @testset "3-bus unbalanced fotr opf with dy transformer" begin
        math = transform_data_model(ravens_ut_trans_2w_dy_lag)
        result = solve_mc_opf(math, FOTRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 467.699; atol=10*200)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), 485.553; atol=10*60)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["3"]["vm"] ./ vbase, [0.92092, 0.91012, 0.90059]; atol=10*3e-1))
    end

    @testset "3-bus unbalanced fotr opf with voltage-dependent loads" begin
        math = transform_data_model(ravens_case3_unbalanced_delta_loads)
        result = solve_mc_opf(math, FOTRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_INFEASIBLE || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 42.0464; atol=10*1)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), 18.1928; atol=10*1)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
   
        @test all(isapprox.(result["bus"]["loadbus"]["vm"] ./ vbase, [0.94105, 0.95942, 0.95876]; atol=10*2e-2))
        @test all(isapprox.(result["bus"]["loadbus"]["va"], [-0.9, -120.3, 120.2]; atol=10*2e-1))
    end
end
