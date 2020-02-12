@info "running current-voltage power flow (pf_iv) tests"

@testset "test pf" begin
    @testset "test opendss pf" begin
        @testset "2-bus diagonal ivr pf" begin
            pmd = PMD.parse_file("../test/data/opendss/case2_diag.dss")
            sol = PMD.run_mc_pf_iv(pmd, PMs.IVRPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
        end

        @testset "3-bus balanced ivr pf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
            sol = PMD.run_mc_pf_iv(pmd, PMs.IVRPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183456; atol=1e-5)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923328; atol=1e-4)
        end
    end
end
