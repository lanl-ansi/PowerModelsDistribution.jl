@info "running current-voltage optimal power flow (opf_iv) tests"

@testset "test current-voltage formulations" begin
    @testset "test IVR opf_iv" begin
        @testset "2-bus diagonal acp opf" begin
            pmd = parse_file("../test/data/opendss/case2_diag.dss")
            sol = run_mc_opf(pmd, IVRPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.018208969542066918; atol = 1e-4)

            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * sol["solution"]["settings"]["sbase"]), 0.018209; atol=1e-5)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * sol["solution"]["settings"]["sbase"]), 0.000208979; atol=1e-5)
        end

        @testset "3-bus balanced acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_balanced.dss")
            sol = run_mc_opf(pmd, IVRPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.018345004773175046; atol = 1e-4)

            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * sol["solution"]["settings"]["sbase"]), 0.018345; atol=1e-6)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * sol["solution"]["settings"]["sbase"]), 0.00919404; atol=1.2e-5)
        end

        @testset "3-bus unbalanced acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_unbalanced.dss")
            sol = run_mc_opf(pmd, IVRPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.021481176584287; atol = 1e-4)


            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * sol["solution"]["settings"]["sbase"]), 0.0214812; atol=1e-6)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * sol["solution"]["settings"]["sbase"]), 0.00927263; atol=1e-5)
        end
    end
end
