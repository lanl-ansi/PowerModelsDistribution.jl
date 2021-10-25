@info "running branch-flow power flow (pf_bf) tests"

@testset "test UBF formulations in pf" begin
    @testset "test lindistflow pf" begin
        @testset "3-bus balanced lpubfdiag pf" begin
            result = solve_mc_pf(case3_balanced, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED

            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 18.0; atol=1e-1)
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), 9.0; atol=1e-1)
        end

        @testset "3-bus unbalanced lpubfdiag pf_bf" begin
            result = solve_mc_pf(case3_unbalanced, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 21.0; atol=1e-1)
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), 9.0; atol=1e-1)
        end
    end

    @testset "test UBF relaxations in pf" begin
        data = deepcopy(case3_unbalanced)
        make_lossless!(data)
        apply_voltage_bounds!(data; vm_lb=0.9, vm_ub=1.1)

        data["settings"]["sbase_default"] = 1.0

        merge!(data["voltage_source"]["source"], Dict{String,Any}(
            "cost_pg_parameters" => [0.0, 1000.0, 0.0],
            "pg_lb" => fill(  0.0, 3),
            "pg_ub" => fill( 10.0, 3),
            "qg_lb" => fill(-10.0, 3),
            "qg_ub" => fill( 10.0, 3),
            )
        )

        for (_,line) in data["line"]
            line["sm_ub"] = fill(10.0, 3)
            delete!(line, "cm_ub")
        end

        @testset "test sdp UBF pf_bf" begin
            @testset "3-bus SDPUBF pf_bf" begin
                result = solve_mc_pf(data, SDPUBFPowerModel, scs_solver)

                @test result["termination_status"] == OPTIMAL
                @test isapprox(result["objective"], 0; atol = 1e-2)
            end
        end

        # @testset "test sdp UBF pf_bf in full matrix form" begin
        #     @testset "3-bus SDPUBFKCLMX pf_bf" begin
        #         result = solve_mc_pf(data, SDPUBFKCLMXPowerModel, scs_solver)
        #
        #         @test result["termination_status"] == OPTIMAL
        #         @test isapprox(result["objective"], 0; atol = 1e-2)
        #     end
        # end


        @testset "test soc UBF pf_bf" begin
            @testset "3-bus SOCNLPUBF pf_bf" begin
                result = solve_mc_pf(data, SOCNLPUBFPowerModel, ipopt_solver)

                @test result["termination_status"] == LOCALLY_SOLVED
                @test isapprox(result["objective"], 0; atol = 1e-1)
            end
            @testset "3-bus SOCConicUBF pf_bf" begin
                result = solve_mc_pf(data, SOCConicUBFPowerModel, scs_solver)

                @test result["termination_status"] == OPTIMAL
                @test isapprox(result["objective"], 0; atol = 1e-2)
            end
        end
    end
end
