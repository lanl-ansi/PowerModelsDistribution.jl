@info "running storage tests"

@testset "test storage opf" begin
    @testset "3-bus balanced battery acp opf" begin
        data = deepcopy(case3_balanced_battery)
        data["settings"]["sbase_default"] = 1e5

        result = solve_mc_opf(case3_balanced_battery, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case3_balanced_battery["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, 0.98697; atol=1e-2))
    end

    @testset "3-bus balanced battery acr opf" begin
        result = solve_mc_opf(case3_balanced_battery, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case3_balanced_battery["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, 0.98697; atol=1e-2))
    end

    @testset "3-bus balanced battery lpubfdiag opf" begin
        result = solve_mc_opf(case3_balanced_battery, LPUBFDiagPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED
        vbase = case3_balanced_battery["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, 0.99767; atol=1e-2))
    end

    @testset "3-bus balanced battery nfa opf" begin
        result = solve_mc_opf(case3_balanced_battery, NFAUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
    end
end

@testset "test storage pf" begin
    @testset "3-bus balanced battery acp pf" begin
        data = deepcopy(case3_balanced_battery)
        data["settings"]["sbase_default"] = 1e5

        result = solve_mc_pf(case3_balanced_battery, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case3_balanced_battery["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, 0.98697; atol=1e-2))
    end

    @testset "3-bus balanced battery acr pf" begin
        result = solve_mc_pf(case3_balanced_battery, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case3_balanced_battery["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, 0.98697; atol=1e-2))
    end
end
