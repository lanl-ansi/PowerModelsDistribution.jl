@info "running branch-flow optimal power flow (opf_bf) tests"

@testset "test distflow formulations in opf" begin
    @testset "test linearised distflow opf_bf" begin
        @testset "3-bus balanced lpubfdiag opf_bf" begin
            result = solve_mc_opf(case3_balanced, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED

            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 18.3456; atol=1)
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), 9.23328; atol=1)
        end

        @testset "3-bus unbalanced lpubfdiag opf_bf" begin
            result = solve_mc_opf(case3_unbalanced, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED

            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 21.4812; atol=1)
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), 9.27263; atol=1)
        end

        @testset "3-bus unbalanced lpubfdiag opf_bf with only two terminals on load bus" begin
            result = solve_mc_opf(case3_unbalanced_missingedge, LPUBFDiagPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test result["termination_status"] == LOCALLY_SOLVED

            vbase = case3_unbalanced_missingedge["settings"]["vbases_default"]["sourcebus"]
            @test all(isapprox.(result["solution"]["bus"]["loadbus"]["vm"] ./ vbase, [0.96038, 0.97866]; atol=1e-4))
        end

        @testset "3-bus unbalanced lpubfdiag opf_bf with delta loads" begin
            data = deepcopy(case3_unbalanced_delta_loads)
            apply_voltage_bounds!(data; vm_lb=0.95, vm_ub=1.05)

            solver = optimizer_with_attributes(
                Ipopt.Optimizer,
                "sb"=>"yes",
                "print_level"=>0,
                "warm_start_init_point"=>"yes",
                "mu_strategy"=>"adaptive"
            )
            result = solve_mc_opf(data, LPUBFDiagPowerModel, solver; solution_processors=[sol_data_model!])

            @test result["termination_status"] == LOCALLY_SOLVED

            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 40.26874; atol=1)
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]),  17.1721; atol=1)

            vbase = case3_unbalanced_delta_loads["settings"]["vbases_default"]["sourcebus"]
            @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, [0.98514,0.98945,0.98929]; atol=5e-2))
            @test all(isapprox.(result["solution"]["bus"]["loadbus"]["vm"] ./ vbase, [0.97007,0.97949,0.97916]; atol=5e-2))
        end

        @testset "3-bus unbalanced fbs opf_bf with switch" begin
            result = solve_mc_opf(case3_unbalanced_switch, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test result["termination_status"] == LOCALLY_SOLVED

            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 21.2194; atol=1)
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), 9.12439; atol=1)

            vbase = case3_unbalanced_switch["settings"]["vbases_default"]["sourcebus"]
            @test all(isapprox.(result["solution"]["bus"]["loadbus"]["vm"] ./ vbase, [0.98102, 0.98922, 0.98692]; atol=9e-2))
            @test all(isapprox.(result["solution"]["bus"]["loadbus"]["va"], [-0.2312, -120.1135, 120.1174]; atol=3e-2))
        end

        @testset "3-bus unbalanced fbs opf_bf with yy transformer" begin
            result = solve_mc_opf(ut_trans_2w_yy, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test result["termination_status"] == LOCALLY_SOLVED

            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 467.547; atol=40)
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), 484.327; atol=20)

            vbase, _ = calc_voltage_bases(ut_trans_2w_yy, ut_trans_2w_yy["settings"]["vbases_default"])
            @test all(isapprox.(result["solution"]["bus"]["3"]["vm"] ./ vbase["3"], [0.87451, 0.8613, 0.85348]; atol=3e-1))
            @test all(isapprox.(result["solution"]["bus"]["3"]["va"], [-0.1, -120.4, 119.8]; atol=2e-1))
        end

        @testset "3-bus unbalanced fbs opf_bf with dy transformer" begin
            result = solve_mc_opf(ut_trans_2w_dy_lag, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test result["termination_status"] == LOCALLY_SOLVED

            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 467.699; atol=200)
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), 485.553; atol=100)

            vbase, _ = calc_voltage_bases(ut_trans_2w_dy_lag, ut_trans_2w_dy_lag["settings"]["vbases_default"])
            @test all(isapprox.(result["solution"]["bus"]["3"]["vm"] ./ vbase["3"], [0.92092, 0.91012, 0.90059]; atol=5e-2))
            @test all(isapprox.(result["solution"]["bus"]["3"]["va"], [-30, -150.4, 89.8]; atol=2e0))
        end

        @testset "3-bus unbalanced fbs opf_bf with voltage-dependent loads" begin
            result = solve_mc_opf(case3_unbalanced_delta_loads, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test result["termination_status"] == LOCALLY_SOLVED

            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 42.0464; atol=1)
            @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), 18.1928; atol=1)

            vbase = case3_unbalanced_delta_loads["settings"]["vbases_default"]["sourcebus"]
            @test all(isapprox.(result["solution"]["bus"]["loadbus"]["vm"] ./ vbase, [0.94105, 0.95942, 0.95876]; atol=2e-3))
            @test all(isapprox.(result["solution"]["bus"]["loadbus"]["va"], [-0.9, -120.3, 120.2]; atol=6e-2))
        end
    end

    @testset "UBF realaxations opf" begin
        data = deepcopy(case3_unbalanced)
        make_lossless!(data; exclude=["line", "linecode"])
        remove_line_limits!(data)
        apply_voltage_bounds!(data; vm_lb=0.9, vm_ub=1.1)

        data["settings"]["sbase_default"] = 1.0

        data["generator"] = Dict{String,Any}(
            "1" => Dict{String,Any}(
                "bus" => "primary",
                "connections" => [1, 2, 3, 4],
                "cost_pg_parameters" => [0.0, 1200.0, 0.0],
                "qg_lb" => fill(0.0, 3),
                "qg_ub" => fill(0.0, 3),
                "pg_ub" => fill(10, 3),
                "pg_lb" => fill(0, 3),
                "configuration" => WYE,
                "status" => ENABLED
            )
        )

        merge!(data["voltage_source"]["source"], Dict{String,Any}(
            "cost_pg_parameters" => [0.0, 1000.0, 0.0],
            "pg_lb" => fill(  0.0, 3),
            "pg_ub" => fill( 10.0, 3),
            "qg_lb" => fill(-10.0, 3),
            "qg_ub" => fill( 10.0, 3),
        ))

        for (_,line) in data["line"]
            line["sm_ub"] = fill(10.0, 3)
        end

        @testset "test sdp distflow opf_bf" begin
            @testset "3-bus SDPUBF opf_bf" begin
                result = solve_mc_opf(data, SDPUBFPowerModel, scs_solver)

                @test result["termination_status"] == OPTIMAL

                @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 21.48; atol = 1e-2)
            end
        end

        # TODO track down why this problem is infeasible in the matrix form (extra Pg,Qg variables?)
        # @testset "test sdp distflow opf_bf in full matrix form" begin
        #     @testset "3-bus SDPUBFKCLMX opf_bf" begin
        #         result = solve_mc_opf(data, SDPUBFKCLMXPowerModel, scs_solver)

        #         @test result["termination_status"] == OPTIMAL
        #         @test isapprox(result["objective"], 21.48; atol = 1e-2)
        #     end
        # end


        @testset "test soc distflow opf_bf" begin
            @testset "3-bus SOCNLPUBF opf_bf" begin
                result = solve_mc_opf(data, SOCNLPUBFPowerModel, ipopt_solver)

                @test result["termination_status"] == LOCALLY_SOLVED

                @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 21.179; atol = 1e-1)
            end
            @testset "3-bus SOCConicUBF opf_bf" begin
                result = solve_mc_opf(data, SOCConicUBFPowerModel, scs_solver)

                @test result["termination_status"] == OPTIMAL || result["termination_status"] == ALMOST_OPTIMAL

                @test isapprox(result["objective"], 21.17; atol = 5e-2)
            end
        end

        @testset "test ubf relaxations with with switches" begin
            data = deepcopy(case3_balanced_switch)
            make_lossless!(data; exclude=["line", "linecode"])
            remove_line_limits!(data)
            apply_voltage_bounds!(data; vm_lb=0.9, vm_ub=1.1)

            data["settings"]["sbase_default"] = 1.0

            data["generator"] = Dict{String,Any}(
                "1" => Dict{String,Any}(
                    "bus" => "primary",
                    "connections" => [1, 2, 3, 4],
                    "cost_pg_parameters" => [0.0, 1200.0, 0.0],
                    "qg_lb" => fill(0.0, 3),
                    "qg_ub" => fill(0.0, 3),
                    "pg_ub" => fill(10, 3),
                    "pg_lb" => fill(0, 3),
                    "configuration" => WYE,
                    "status" => ENABLED
                )
            )

            merge!(data["voltage_source"]["source"], Dict{String,Any}(
                "cost_pg_parameters" => [0.0, 1000.0, 0.0],
                "pg_lb" => fill(  0.0, 3),
                "pg_ub" => fill( 10.0, 3),
                "qg_lb" => fill(-10.0, 3),
                "qg_ub" => fill( 10.0, 3),
            ))

            for (_,line) in data["line"]
                line["sm_ub"] = fill(10.0, 3)
            end

            @testset "test SOCNLPUBF opf with switches" begin
                result = solve_mc_opf(data, SOCNLPUBFPowerModel, ipopt_solver)

                @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
                @test all(isapprox.(result["solution"]["switch"]["ohline"]["pf"], [6.0, 6.0, 6.0]; atol=1e-1))
                @test isapprox(result["objective"], 18.1824; atol=2e-1)
            end

            @testset "test SOCConicUBF opf with switches" begin
                result = solve_mc_opf(data, SOCConicUBFPowerModel, scs_solver)

                @test result["termination_status"] == OPTIMAL || result["termination_status"] == ALMOST_OPTIMAL
                @test all(isapprox.(result["solution"]["switch"]["ohline"]["pf"], [6.0, 6.0, 6.0]; atol=1e-1))
                @test isapprox(result["objective"], 18.1824; atol=2e-1)
            end
        end
    end
end
