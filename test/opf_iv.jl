@info "running current-voltage optimal power flow (opf_iv) tests"

@testset "test current-voltage formulations" begin
    @testset "test IVR opf_iv" begin
        @testset "2-bus diagonal ivr opf" begin
            result = solve_mc_opf(case2_diag, IVRUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED

            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 18.209; atol=1e-2)
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]),  0.208; atol=1e-2)
        end

        @testset "3-bus balanced ivr opf" begin
            result = solve_mc_opf(case3_balanced, IVRUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 18.345; atol=1e-2)
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]),  9.194; atol=1e-2)
        end

        @testset "3-bus unbalanced ivr opf" begin
            result = solve_mc_opf(case3_unbalanced, IVRUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED

            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 21.4812; atol=1e-2)
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), 9.27263; atol=1e-2)
        end

        @testset "ivr opf power variable expressions" begin
            math = transform_data_model(IEEE13_Assets)

            pm = instantiate_mc_model(IEEE13_Assets, IVRUPowerModel, build_mc_opf)

            @test length(var(pm, :p)) == length(math["branch"])*2 && length(var(pm, :q)) == length(math["branch"])*2
            @test length(var(pm, :pt)) == length(math["transformer"])*2 && length(var(pm, :qt)) == length(math["transformer"])*2
            @test length(var(pm, :psw)) == length(math["switch"])*2 && length(var(pm, :qsw)) == length(math["switch"])*2
        end
    end
end
