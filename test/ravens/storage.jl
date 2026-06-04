@info "running storage tests"

@testset "test storage opf" begin
    @testset "3-bus balanced battery acp opf" begin
        math = transform_data_model(ravens_case3_balanced_battery)
        result = solve_mc_opf(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["primary"]["vm"] ./ vbase, 0.991111; atol=1e-2))
        @test isapprox(sum(result["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end

    @testset "3-bus balanced battery acp opf - time_elapsed::Int" begin
        math = transform_data_model(ravens_case3_balanced_battery)
        case = deepcopy(math)
        case["time_elapsed"] = 1

        result = solve_mc_opf(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["primary"]["vm"] ./ vbase, 0.991111; atol=1e-2))
        @test isapprox(sum(result["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end

    @testset "3-bus balanced battery acr opf" begin
        math = transform_data_model(ravens_case3_balanced_battery)
        result = solve_mc_opf(math, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["primary"]["vm"] ./ vbase, 0.991111; atol=1e-2))
        @test isapprox(sum(result["storage"]["s1"]["ps"]), 0.0; atol=1e-4)#changed to zero
    end

    @testset "3-bus balanced battery lpubfdiag opf" begin
        math = transform_data_model(ravens_case3_balanced_battery)
        result = solve_mc_opf(math, LPUBFDiagPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["primary"]["vm"] ./ vbase, 0.991111; atol=1e-2))
        @test isapprox(sum(result["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end

    @testset "3-bus balanced battery nfa opf" begin
        math = transform_data_model(ravens_case3_balanced_battery)
        result = solve_mc_opf(math, NFAUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)
        @test isapprox(sum(result["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end
end

@testset "test storage pf" begin
    @testset "3-bus balanced battery acp pf" begin
        math = transform_data_model(ravens_case3_balanced_battery)
        result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["primary"]["vm"] ./ vbase, 0.991111; atol=1e-2))
        @test all(isapprox.(result["bus"]["primary"]["va"], [0.03, -119.97, 120.03]; atol=1e-2))
        @test isapprox(sum(result["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end

    @testset "3-bus balanced battery acr pf" begin
        math = transform_data_model(ravens_case3_balanced_battery)
        result = solve_mc_pf(math, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        # Test is numerically unstable (fails on only some OSes and some versions of Julia)
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["primary"]["vm"] ./ vbase, 0.991111; atol=1e-2))
        @test all(isapprox.(result["bus"]["primary"]["va"], [0.03, -119.97, 120.03]; atol=1e-2))
        @test isapprox(sum(result["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end
end
