@info "running branch-flow power flow (pf_bf) tests"

@testset "test distflow formulations in pf" begin
    @testset "test linearised distflow pf_bf" begin
        @testset "5-bus lpubfdiag opf_bf" begin
            mp_data = PowerModels.parse_file("../test/data/matpower/case5.m")
            make_multiconductor!(mp_data, 3)
            result = solve_mc_pf(mp_data, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0; atol = 1e0)
            # @test isapprox(result["solution"]["bus"]["3"]["vm"], 0.911466*[1,1,1]; atol = 1e-3)
            vm = calc_vm_w(result, "3")
            @test isapprox(vm, [1,1,1]; atol = 1e-3)

        end

        @testset "3-bus balanced lpubfdiag pf_bf" begin
            pmd = parse_file("../test/data/opendss/case3_balanced.dss"; data_model=MATHEMATICAL)
            sol = solve_mc_pf(pmd, LPUBFDiagPowerModel, ipopt_solver)

            @test sol["termination_status"] == LOCALLY_SOLVED
            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * baseMVA), 0.0183456; atol=2e-3)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * baseMVA), 0.00923328; atol=2e-3)
        end

        @testset "3-bus unbalanced lpubfdiag pf_bf" begin
            pmd = parse_file("../test/data/opendss/case3_unbalanced.dss"; data_model=MATHEMATICAL)
            sol = solve_mc_pf(pmd, LPUBFDiagPowerModel, ipopt_solver)

            @test sol["termination_status"] == LOCALLY_SOLVED
            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * baseMVA), 0.0214812; atol=2e-3)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * baseMVA), 0.00927263; atol=2e-3)
        end
    end

    @testset "test linearised distflow pf_bf in diagonal matrix form" begin
        @testset "5-bus lpdiagubf pf_bf" begin
            mp_data = PowerModels.parse_file("../test/data/matpower/case5.m")
            make_multiconductor!(mp_data, 3)
            result = solve_mc_pf(mp_data, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0; atol = 1e0)
        end
    end

    @testset "test linearised distflow pf_bf in full matrix form" begin
        @testset "5-bus lpfullubf pf_bf" begin
            mp_data = PowerModels.parse_file("../test/data/matpower/case5.m")
            make_multiconductor!(mp_data, 3)
            result = solve_mc_pf(mp_data, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0; atol = 1e0)
        end
    end

    data = parse_file("../test/data/opendss/case3_unbalanced.dss"; transformations=[make_lossless!])
    data["settings"]["sbase_default"] = 0.001 * 1e3
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
    end

    data = transform_data_model(data)

    for (_,bus) in data["bus"]
        if bus["name"] != "sourcebus"
            bus["vmin"] = fill(0.9, 3)
            bus["vmax"] = fill(1.1, 3)
        end
    end

    @testset "test sdp distflow pf_bf" begin
        @testset "3-bus SDPUBF pf_bf" begin
            result = solve_mc_pf(data, SDPUBFPowerModel, scs_solver)

            @test result["termination_status"] == OPTIMAL
            @test isapprox(result["objective"], 0; atol = 1e-2)
        end
    end

    # @testset "test sdp distflow pf_bf in full matrix form" begin
    #     @testset "3-bus SDPUBFKCLMX pf_bf" begin
    #         result = solve_mc_pf(data, SDPUBFKCLMXPowerModel, scs_solver)
    #
    #         @test result["termination_status"] == OPTIMAL
    #         @test isapprox(result["objective"], 0; atol = 1e-2)
    #     end
    # end


    @testset "test soc distflow pf_bf" begin
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
