@info "running storage tests"

@testset "test storage opf" begin
    @testset "3-bus balanced battery acp opf" begin
        result = solve_mc_opf(case3_balanced_battery, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case3_balanced_battery["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, 0.991111; atol=1e-2))
        @test isapprox(sum(result["solution"]["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end

    @testset "3-bus balanced battery acp opf - time_elapsed::Int" begin
        case = deepcopy(case3_balanced_battery)
        case["time_elapsed"] = 1

        result = solve_mc_opf(case3_balanced_battery, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case3_balanced_battery["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, 0.991111; atol=1e-2))
        @test isapprox(sum(result["solution"]["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end

    @testset "3-bus balanced battery acr opf" begin
        result = solve_mc_opf(case3_balanced_battery, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case3_balanced_battery["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, 0.991111; atol=1e-2))
        @test isapprox(sum(result["solution"]["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end

    @testset "3-bus balanced battery lpubfdiag opf" begin
        result = solve_mc_opf(case3_balanced_battery, LPUBFDiagPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED
        vbase = case3_balanced_battery["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, 0.991111; atol=1e-2))
        @test isapprox(sum(result["solution"]["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end

    @testset "3-bus balanced battery nfa opf" begin
        result = solve_mc_opf(case3_balanced_battery, NFAUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(sum(result["solution"]["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end
end

@testset "test storage pf" begin
    @testset "3-bus balanced battery acp pf" begin
        result = solve_mc_pf(case3_balanced_battery, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case3_balanced_battery["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, 0.991111; atol=1e-2))
        @test all(isapprox.(result["solution"]["bus"]["primary"]["va"], [0.03, -119.97, 120.03]; atol=1e-2))
        @test isapprox(sum(result["solution"]["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end

    @testset "3-bus balanced battery acr pf" begin
        result = solve_mc_pf(case3_balanced_battery, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        # Test is numerically unstable (fails on only some OSes and some versions of Julia)
        vbase = case3_balanced_battery["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, 0.991111; atol=1e-2))
        @test all(isapprox.(result["solution"]["bus"]["primary"]["va"], [0.03, -119.97, 120.03]; atol=1e-2))
        @test isapprox(sum(result["solution"]["storage"]["s1"]["ps"]), -5.0; atol=1e-4)
    end
end
