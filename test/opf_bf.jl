@info "running branch-flow optimal power flow (opf_bf) tests"

@testset "test distflow formulations" begin
    @testset "test linearised distflow opf_bf" begin
        @testset "5-bus lpubfdiag opf_bf" begin
            mp_data = PowerModels.parse_file("../test/data/matpower/case5.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf_bf(mp_data, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 44880; atol = 1e0)
            # @test isapprox(result["solution"]["bus"]["3"]["vm"], 0.911466*[1,1,1]; atol = 1e-3)
            vm = calc_vm_w(result, "3")
            @test isapprox(vm, 0.911466*[1,1,1]; atol = 1e-3)

        end

        @testset "3-bus balanced lpubfdiag opf_bf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
            sol = PMD.run_mc_opf_bf(pmd, LPUBFDiagPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183456; atol=2e-3)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923328; atol=2e-3)
        end

        @testset "3-bus unbalanced lpubfdiag opf_bf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced.dss")
            sol = PMD.run_mc_opf_bf(pmd, LPUBFDiagPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214812; atol=2e-3)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00927263; atol=2e-3)
        end
    end

    @testset "test linearised distflow opf_bf in diagonal matrix form" begin
        @testset "5-bus lpdiagubf opf_bf" begin
            mp_data = PowerModels.parse_file("../test/data/matpower/case5.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf_bf(mp_data, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 44880; atol = 1e0)
        end
    end

    @testset "test linearised distflow opf_bf in full matrix form" begin
        @testset "5-bus lpfullubf opf_bf" begin
            mp_data = PowerModels.parse_file("../test/data/matpower/case5.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf_bf(mp_data, LPUBFDiagPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 44880; atol = 1e0)
        end
    end

    function case3_unbalanced_with_bounds()
        file = "../test/data/opendss/case3_unbalanced.dss"
        data = PowerModelsDistribution.parse_file(file)

        baseMVA = data["baseMVA"]
        #set some reasonable bounds
        data["gen"]["1"]["cost"] =  100*data["gen"]["1"]["cost"]
        data["gen"]["1"]["pmin"] = 0*[1.0, 1.0, 1.0]/baseMVA
        data["gen"]["1"]["pmax"] = [1.0, 1.0, 1.0]/baseMVA
        data["gen"]["1"]["qmin"] = -1*[1.0, 1.0, 1.0]/baseMVA
        data["gen"]["1"]["qmax"] =  [1.0, 1.0, 1.0]/baseMVA

        #improve numerical properties by removing source bus impedance
        data["bus"]["1"]["bus_type"] = 3
        data["bus"]["1"]["vm"] = data["bus"]["4"]["vm"]
        data["bus"]["1"]["vmin"] = data["bus"]["4"]["vmin"]
        data["bus"]["1"]["vmax"] = data["bus"]["4"]["vmax"]
        delete!(data["branch"], "3")
        delete!(data["bus"], "4")
        data["gen"]["1"]["gen_bus"] = 1

        for (n, branch) in (data["branch"])
                branch["rate_a"] = [1.0, 1.0, 1.0]/baseMVA
        end
        return data
    end

    @testset "test sdp distflow opf_bf" begin
        @testset "3-bus SDPUBF opf_bf" begin
            data = case3_unbalanced_with_bounds()
            result = run_mc_opf_bf(data, SDPUBFPowerModel, scs_solver)

            @test result["termination_status"] == PMs.OPTIMAL
            @test isapprox(result["objective"], 2.148; atol = 1e-3)
        end
    end

    @testset "test sdp distflow opf_bf in full matrix form" begin
        @testset "3-bus SDPUBFKCLMX opf_bf" begin
            data = case3_unbalanced_with_bounds()
            result = run_mc_opf_bf(data, SDPUBFKCLMXPowerModel, scs_solver)

            @test result["termination_status"] == PMs.OPTIMAL
            @test isapprox(result["objective"], 2.148; atol = 1e-3)
        end
    end


    @testset "test soc distflow opf_bf" begin
        @testset "3-bus SOCNLPUBF opf_bf" begin
            data = case3_unbalanced_with_bounds()
            result = run_mc_opf_bf(data, SOCNLPUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 2.117; atol = 1e-2)
        end
        @testset "3-bus SOCConicUBF opf_bf" begin
            data = case3_unbalanced_with_bounds()
            result = run_mc_opf_bf(data, SOCConicUBFPowerModel, scs_solver)

            @test result["termination_status"] == PMs.ALMOST_OPTIMAL
            @test isapprox(result["objective"], 2.117; atol = 1e-3)
        end
    end
end
