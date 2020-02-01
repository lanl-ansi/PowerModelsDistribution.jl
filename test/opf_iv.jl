@info "running current-voltage optimal power flow (opf_iv) tests"

@testset "test current-voltage formulations" begin
    @testset "test IVR opf_iv" begin
        @testset "2-bus diagonal acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case2_diag.dss")
            sol = PMD.run_mc_opf_iv(pmd, PMs.IVRPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.018208969542066918; atol = 1e-4)

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
        end

        @testset "3-bus balanced acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
            sol = PMD.run_mc_opf_iv(pmd, PMs.IVRPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.018345004773175046; atol = 1e-4)

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018345; atol=1e-6)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00919404; atol=1.2e-5)
        end

        @testset "3-bus unbalanced acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced.dss")
            sol = PMD.run_mc_opf_iv(pmd, PMs.IVRPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.021481176584287; atol = 1e-4)


            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214812; atol=1e-6)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00927263; atol=1e-5)
        end
    end
end
