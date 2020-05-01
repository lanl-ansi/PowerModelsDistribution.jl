@info "running branch-flow optimal power flow (opf_bf) tests"

@testset "test distflow formulations" begin
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
end
