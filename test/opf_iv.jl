@info "running current-voltage optimal power flow (opf_iv) tests"

@testset "test current-voltage formulations" begin
    @testset "test IVR opf_iv" begin
        @testset "2-bus diagonal acp opf" begin
            pmd = parse_file("../test/data/opendss/case2_diag.dss")
            sol = solve_mc_opf(pmd, IVRUPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.018208969542066918; atol = 1e-4)

            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA), 0.018209; atol=1e-5)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA), 0.000208979; atol=1e-5)
        end

        @testset "3-bus balanced acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_balanced.dss")
            sol = solve_mc_opf(pmd, IVRUPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.018345004773175046; atol = 1e-4)

            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA), 0.018345; atol=1e-6)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA), 0.00919404; atol=1.2e-5)
        end

        @testset "3-bus unbalanced acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_unbalanced.dss")
            sol = solve_mc_opf(pmd, IVRUPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.021481176584287; atol = 1e-4)

            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA), 0.0214812; atol=1e-6)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA), 0.00927263; atol=1e-5)
        end

        @testset "5-bus storage matpower mn ivr opf" begin
            case5 = PM.parse_file("../test/data/matpower/case5.m")
            make_multiconductor!(case5, 3)
            case5_mn = InfrastructureModels.replicate(case5, 3, Set(["per_unit"]))
            result = solve_mn_mc_opf(case5_mn, IVRUPowerModel, ipopt_solver)
            @test result["termination_status"] == LOCALLY_SOLVED
        end
    end
end
