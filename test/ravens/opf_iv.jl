@info "running current-voltage optimal power flow (opf_iv) tests"

@testset "test current-voltage formulations" begin
    @testset "test IVR opf_iv" begin
        @testset "2-bus diagonal ivr opf" begin
            math = transform_data_model(ravens_case2_diag)
            result = solve_mc_opf(math, IVRUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

            result = transform_solution(result["solution"], math, make_si=true)

            @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 18.209; atol=1e-2)
            @test isapprox(sum(result["voltage_source"]["source"]["qg"]),  0.208; atol=1e-2)
        end

        @testset "3-bus balanced ivr opf" begin
            math = transform_data_model(ravens_case3_balanced)
            result = solve_mc_opf(math, IVRUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

            result = transform_solution(result["solution"], math, make_si=true)

            @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 18.345; atol=1e-2)
            @test isapprox(sum(result["voltage_source"]["source"]["qg"]),  9.194; atol=1e-2)
        end

        @testset "3-bus unbalanced ivr opf" begin
            math = transform_data_model(ravens_case3_unbalanced)
            result = solve_mc_opf(math, IVRUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

            result = transform_solution(result["solution"], math, make_si=true)

            @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 21.4812; atol=1e-2)
            @test isapprox(sum(result["voltage_source"]["source"]["qg"]), 9.27263; atol=1e-2)
        end

        @testset "ivr opf power variable expressions" begin
            math = transform_data_model(ravens_IEEE13_Assets)

            pm = instantiate_mc_model(IEEE13_Assets, IVRUPowerModel, build_mc_opf)
            
            @test length(var(pm, :p)) == (length(math["branch"])+1)*2 && length(var(pm, :q)) == (length(math["branch"])+1)*2
            @test length(var(pm, :pt)) == length(math["transformer"])*2 && length(var(pm, :qt)) == length(math["transformer"])*2
            @test length(var(pm, :psw)) == length(math["switch"])*2 && length(var(pm, :qsw)) == length(math["switch"])*2
        end
    end
end

