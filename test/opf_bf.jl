@info "running branch-flow optimal power flow (opf_bf) tests"

@testset "test distflow formulations in opf" begin
    case5 = PowerModels.parse_file("../test/data/matpower/case5.m")
    make_multiconductor!(case5, 3)

    @testset "test linearised distflow opf_bf" begin
        @testset "5-bus lpubfdiag opf_bf" begin
            result = run_mc_opf(case5, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 44880; atol = 1e0)
            # @test isapprox(result["solution"]["bus"]["3"]["vm"], 0.911466*[1,1,1]; atol = 1e-3)
            vm = calc_vm_w(result, "3")
            @test isapprox(vm, 0.911466*[1,1,1]; atol = 1e-3)

        end

        @testset "3-bus balanced lpubfdiag opf_bf" begin
            pmd = parse_file("../test/data/opendss/case3_balanced.dss")
            sol = run_mc_opf(pmd, LPUBFDiagPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * sol["solution"]["settings"]["sbase"]), 0.0183456; atol=2e-3)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * sol["solution"]["settings"]["sbase"]), 0.00923328; atol=2e-3)
        end

        @testset "3-bus unbalanced lpubfdiag opf_bf" begin
            pmd = parse_file("../test/data/opendss/case3_unbalanced.dss")
            sol = run_mc_opf(pmd, LPUBFDiagPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * sol["solution"]["settings"]["sbase"]), 0.0214812; atol=2e-3)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * sol["solution"]["settings"]["sbase"]), 0.00927263; atol=2e-3)
        end
    end

    @testset "test linearised distflow opf_bf in diagonal matrix form" begin
        @testset "5-bus lpdiagubf opf_bf" begin
            result = run_mc_opf(case5, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 44880; atol = 1e0)
        end
    end

    @testset "test linearised distflow opf_bf in full matrix form" begin
        @testset "5-bus lpfullubf opf_bf" begin
            result = run_mc_opf(case5, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 44880; atol = 1e0)
        end
    end

    data = parse_file("../test/data/opendss/case3_unbalanced.dss"; transformations=[make_lossless!])
    data["settings"]["sbase_default"] = 0.001 * 1e3
    data["generator"] = Dict{Any,Any}(
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
        "qg_ub" => fill( 10.0, 3)
    ))

    for (_,line) in data["line"]
        line["sm_ub"] = fill(10.0, 3)
    end

    data = transform_data_model(data)

    for (_,bus) in data["bus"]
        if bus["name"] != "sourcebus"
            bus["vmin"] = fill(0.9, 3)
            bus["vmax"] = fill(1.1, 3)
            bus["vm"] = fill(1.0, 3)
            bus["va"] = deg2rad.([0., -120, 120])
        end
    end

    @testset "test sdp distflow opf_bf" begin
        @testset "3-bus SDPUBF opf_bf" begin
            result = run_mc_opf(data, SDPUBFPowerModel, scs_solver)

            @test result["termination_status"] == OPTIMAL
            @test isapprox(result["objective"], 21.48; atol = 1e-2)
        end
    end

    @testset "test sdp distflow opf_bf in full matrix form" begin
        @testset "3-bus SDPUBFKCLMX opf_bf" begin
            result = run_mc_opf(data, SDPUBFKCLMXPowerModel, scs_solver)

            @test result["termination_status"] == OPTIMAL
            @test isapprox(result["objective"], 21.48; atol = 1e-2)
        end
    end


    @testset "test soc distflow opf_bf" begin
        @testset "3-bus SOCNLPUBF opf_bf" begin
            result = run_mc_opf(data, SOCNLPUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 21.179; atol = 1e-1)
        end
        # @testset "3-bus SOCConicUBF opf_bf" begin
        #     result = run_mc_opf(data, SOCConicUBFPowerModel, scs_solver)
        #
        #     @test result["termination_status"] == ALMOST_OPTIMAL
        #     @test isapprox(result["objective"], 21.17; atol = 1e-2)
        # end
    end
end
