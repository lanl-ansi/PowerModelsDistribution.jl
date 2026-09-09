@info "running branch-flow optimal power flow (opf_bf) tests"

@testset "test distflow formulations in opf" begin
    @testset "test linearised distflow opf_bf" begin
        @testset "3-bus balanced lpubfdiag opf_bf" begin
            math = transform_data_model(ravens_case3_balanced)
            result = solve_mc_opf(math, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

            
            result = transform_solution(result["solution"], math, make_si=true)

            @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 18.3456; atol=1)
            @test isapprox(sum(result["voltage_source"]["source"]["qg"]), 9.23328; atol=1)
        end

        @testset "3-bus unbalanced lpubfdiag opf_bf" begin
            math = transform_data_model(ravens_case3_unbalanced)
            result = solve_mc_opf(math, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

            
            result = transform_solution(result["solution"], math, make_si=true)

            @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 21.4812; atol=1)
            @test isapprox(sum(result["voltage_source"]["source"]["qg"]), 9.27263; atol=1)
        end

        @testset "3-bus unbalanced lpubfdiag opf_bf with only two terminals on load bus" begin
            math = transform_data_model(ravens_case3_unbalanced_missingedge)
            result = solve_mc_opf(math, LPUBFDiagPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

            
            result = transform_solution(result["solution"], math, make_si=true)

            vbase = case3_unbalanced_missingedge["settings"]["vbases_default"]["sourcebus"]
            @test all(isapprox.(result["bus"]["loadbus"]["vm"] ./ vbase, [0.96038, 0.97866]; atol=1e-4))
        end

        # @testset "3-bus unbalanced lpubfdiag opf_bf with delta loads" begin
        #     data = transform_data_model(ravens_case3_unbalanced_delta_loads)
        #     apply_voltage_bounds_math!(data; vm_lb=0.95, vm_ub=1.05)

        #     # to fix unbounded objective ipopt issue
        #     println("<debug> "*string(data))

        #     data["voltage_source"] = Dict{String, Dict{String, Any}}()
        #     data["voltage_source"]["source"] = Dict{String, Any}()
        #     data["voltage_source"]["source"]["pg_lb"] = zeros(3)
        #     data["voltage_source"]["source"]["qg_lb"] = zeros(3)
        #     data["voltage_source"]["source"]["pg_ub"] = fill(50.0, 3)

        #     result = solve_mc_opf(data, LPUBFDiagPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        #     @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

            
        #     result = transform_solution(result["solution"], math, make_si=true)

        #     @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 40.26874; atol=1)
        #     @test isapprox(sum(result["voltage_source"]["source"]["qg"]),  17.1721; atol=1)

            
        #     sourcebus_id = math["bus_lookup"]["sourcebus"]
        #     vbase = math["bus"][string(sourcebus_id)]["vbase"]

        #     @test all(isapprox.(result["bus"]["primary"]["vm"] ./ vbase, [0.98514,0.98945,0.98929]; atol=5e-2))
        #     @test all(isapprox.(result["bus"]["loadbus"]["vm"] ./ vbase, [0.97007,0.97949,0.97916]; atol=5e-2))
        # end

        @testset "3-bus unbalanced fbs opf_bf with switch" begin
            math = transform_data_model(ravens_case3_unbalanced_switch)
            result = solve_mc_opf(math, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

            
            result = transform_solution(result["solution"], math, make_si=true)

            @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 21.2194; atol=1)
            @test isapprox(sum(result["voltage_source"]["source"]["qg"]), 9.12439; atol=1)

            
            sourcebus_id = math["bus_lookup"]["sourcebus"]
            vbase = math["bus"][string(sourcebus_id)]["vbase"]

            @test all(isapprox.(result["bus"]["loadbus"]["vm"] ./ vbase, [0.98102, 0.98922, 0.98692]; atol=9e-2))
            @test all(isapprox.(result["bus"]["loadbus"]["va"], [-0.2312, -120.1135, 120.1174]; atol=3e-2))
        end

        @testset "3-bus unbalanced fbs opf_bf with yy transformer" begin
            math = transform_data_model(ravens_ut_trans_2w_yy)
            result = solve_mc_opf(math, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test_broken result["termination_status"] == LOCALLY_SOLVED

            
            result = transform_solution(result["solution"], math, make_si=true)

            @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 467.547; atol=200)
            @test isapprox(sum(result["voltage_source"]["source"]["qg"]), 484.327; atol=150)

            sourcebus_id = math["bus_lookup"]["sourcebus"]
            vbase = math["bus"][string(sourcebus_id)]["vbase"]
            
            @test all(isapprox.(result["bus"]["3"]["vm"] ./ vbase, [0.4366, 0.446, 0.46457]; atol=3e-1))
            @test_broken all(isapprox.(result["bus"]["3"]["va"], [-0.1, -120.4, 119.8]; atol=2))
        end

        @testset "3-bus unbalanced fbs opf_bf with dy transformer" begin
            math = transform_data_model(ravens_ut_trans_2w_dy_lag)
            result = solve_mc_opf(math, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test_broken result["termination_status"] == LOCALLY_SOLVED

            
            result = transform_solution(result["solution"], math, make_si=true)

            @test_broken isapprox(sum(result["voltage_source"]["source"]["pg"]), 467.699; atol=200)
            @test_broken isapprox(sum(result["voltage_source"]["source"]["qg"]), 485.553; atol=150)

            sourcebus_id = math["bus_lookup"]["sourcebus"]
            vbase = math["bus"][string(sourcebus_id)]["vbase"]
            
            @test_broken all(isapprox.(result["bus"]["3"]["vm"] ./ vbase, [0.3003, 0.3002, 0.3002]; atol=5e-2))
            println(result["bus"]["3"]["va"]-[-30, -150.4, 89.8])
            println(result["bus"]["3"]["va"])
            @test_broken all(isapprox.(result["bus"]["3"]["va"], [-30, -150.4, 89.8]; atol=2))
        end

        @testset "3-bus unbalanced fbs opf_bf with voltage-dependent loads" begin
            math = transform_data_model(ravens_case3_unbalanced_delta_loads)
            result = solve_mc_opf(math, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test_broken result["termination_status"] == LOCALLY_SOLVED #LOCALLY_INFEASIBLE

            
            result = transform_solution(result["solution"], math, make_si=true)

            @test_broken isapprox(sum(result["voltage_source"]["source"]["pg"]), 42.0464; atol=1) #21.7 == 42
            @test_broken isapprox(sum(result["voltage_source"]["source"]["qg"]), 18.1928; atol=1) #11.6 = 18

            sourcebus_id = math["bus_lookup"]["sourcebus"]
            vbase = math["bus"][string(sourcebus_id)]["vbase"]
            @test all(isapprox.(result["bus"]["loadbus"]["vm"] ./ vbase, [0.9512, 0.9964, 0.9936]; atol=9e-1))
            @test_broken all(isapprox.(result["bus"]["loadbus"]["va"], [-0.3733, -120.22, 120.06]; atol=6e-2)) #all(isapprox.(((result["bus"])["loadbus"])["va"], [-0.3733, -120.22, 120.06]; atol = 0.06))
        end
    end

    # @testset "UBF realaxations opf" begin
    #     dara = transform_data_model(case3_unbalanced)
    #     make_lossless!(data; exclude=["line", "linecode"])
    #     remove_line_limits!(data)
    #     apply_voltage_bounds!(data; vm_lb=0.9, vm_ub=1.1)

    #     data["settings"]["sbase_default"] = 1.0

    #     data["generator"] = Dict{String,Any}(
    #         "1" => Dict{String,Any}(
    #             "bus" => "primary",
    #             "connections" => [1, 2, 3, 4],
    #             "cost_pg_parameters" => [0.0, 1200.0, 0.0],
    #             "qg_lb" => fill(0.0, 3),
    #             "qg_ub" => fill(0.0, 3),
    #             "pg_ub" => fill(10, 3),
    #             "pg_lb" => fill(0, 3),
    #             "configuration" => WYE,
    #             "status" => ENABLED
    #         )
    #     )

    #     merge!(data["voltage_source"]["source"], Dict{String,Any}(
    #         "cost_pg_parameters" => [0.0, 1000.0, 0.0],
    #         "pg_lb" => fill(  0.0, 3),
    #         "pg_ub" => fill( 10.0, 3),
    #         "qg_lb" => fill(-10.0, 3),
    #         "qg_ub" => fill( 10.0, 3),
    #     ))

    #     for (_,line) in data["line"]
    #         line["sm_ub"] = fill(10.0, 3)
    #     end

    #     @testset "test sdp distflow opf_bf" begin
    #         @testset "3-bus SDPUBF opf_bf" begin
    #             result = solve_mc_opf(data, SDPUBFPowerModel, scs_solver)

    #             @test result["termination_status"] == OPTIMAL

    #             @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 21.48; atol = 1e-2)
    #         end
    #     end

    #     # TODO track down why this problem is infeasible in the matrix form (extra Pg,Qg variables?)
    #     # @testset "test sdp distflow opf_bf in full matrix form" begin
    #     #     @testset "3-bus SDPUBFKCLMX opf_bf" begin
    #     #         result = solve_mc_opf(data, SDPUBFKCLMXPowerModel, scs_solver)

    #     #         @test result["termination_status"] == OPTIMAL
    #     #         @test isapprox(result["objective"], 21.48; atol = 1e-2)
    #     #     end
    #     # end


    #     @testset "test soc distflow opf_bf" begin
    #         @testset "3-bus SOCNLPUBF opf_bf" begin
    #             result = solve_mc_opf(data, SOCNLPUBFPowerModel, ipopt_solver)

    #             @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

    #             @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 21.179; atol = 1e-1)
    #         end
    #         @testset "3-bus SOCConicUBF opf_bf" begin
    #             result = solve_mc_opf(data, SOCConicUBFPowerModel, scs_solver)

    #             @test_skip result["termination_status"] == OPTIMAL || result["termination_status"] == ALMOST_OPTIMAL

    #             @test isapprox(result["objective"], 21.17; atol = 5e-2)
    #         end
    #     end

        # @testset "test ubf relaxations with with switches" begin
        #     data = transform_data_model(ravens_case4_phase_drop)
        #     make_lossless!(data; exclude=["line", "linecode"])
        #     remove_line_limits!(data)
        #     apply_voltage_bounds!(data; vm_lb=0.9, vm_ub=1.1)

        #     data["settings"]["sbase_default"] = 1.0

        #     data["generator"] = Dict{String,Any}(
        #         "1" => Dict{String,Any}(
        #             "bus" => "primary",
        #             "connections" => [1, 2, 3, 4],
        #             "cost_pg_parameters" => [0.0, 1200.0, 0.0],
        #             "qg_lb" => fill(0.0, 3),
        #             "qg_ub" => fill(0.0, 3),
        #             "pg_ub" => fill(10, 3),
        #             "pg_lb" => fill(0, 3),
        #             "configuration" => WYE,
        #             "status" => ENABLED
        #         )
        #     )

        #     merge!(data["voltage_source"]["source"], Dict{String,Any}(
        #         "cost_pg_parameters" => [0.0, 1000.0, 0.0],
        #         "pg_lb" => fill(  0.0, 3),
        #         "pg_ub" => fill( 10.0, 3),
        #         "qg_lb" => fill(-10.0, 3),
        #         "qg_ub" => fill( 10.0, 3),
        #     ))

        #     for (_,line) in data["line"]
        #         line["sm_ub"] = fill(10.0, 3)
        #     end

        #     @testset "test SOCNLPUBF opf with switches" begin
        #         result = solve_mc_opf(data, SOCNLPUBFPowerModel, ipopt_solver)

        #         @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
        #         @test all(isapprox.(result["switch"]["ohline"]["pf"], [6.0, 6.0, 6.0]; atol=1e-1))
        #         @test isapprox(result["objective"], 18.1824; atol=2e-1)
        #     end

        #     @testset "test SOCConicUBF opf with switches" begin
        #         result = solve_mc_opf(data, SOCConicUBFPowerModel, scs_solver)

        #         @test result["termination_status"] == OPTIMAL || result["termination_status"] == ALMOST_OPTIMAL
        #         @test all(isapprox.(result["switch"]["ohline"]["pf"], [6.0, 6.0, 6.0]; atol=1e-1))
        #         @test isapprox(result["objective"], 18.1824; atol=2e-1)
        #     end
        # end
    # end
end
