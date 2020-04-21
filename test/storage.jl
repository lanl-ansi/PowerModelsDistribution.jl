@info "running storage optimal power flow tests"

@testset "test storage opf" begin
    mp_data = PM.parse_file("../test/data/matpower/case5_strg.m")
    make_multiconductor!(mp_data, 3)

    @testset "5-bus storage acp opf_strg" begin
        result = run_mc_opf(mp_data, ACPPowerModel, ipopt_solver; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 52299.2; atol = 1e0)

        @test isapprox(result["solution"]["storage"]["1"]["se"],  0.0; atol=1e0)
        @test isapprox(result["solution"]["storage"]["2"]["se"],  0.0; atol=1e0)

        @test all(isapprox.(result["solution"]["storage"]["1"]["ps"], -0.0596928; atol=1e-3))
        @test all(isapprox.(result["solution"]["storage"]["2"]["ps"], -0.0794600; atol=1e-3))
    end

    @testset "5-bus storage dcp opf_strg" begin
        result = run_mc_opf(mp_data, DCPPowerModel, ipopt_solver; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 52059.6; atol = 1e0)

        @test isapprox(result["solution"]["storage"]["1"]["se"],  0.0; atol = 1e0)
        @test isapprox(result["solution"]["storage"]["2"]["se"],  0.0; atol = 1e0)

        @test all(isapprox.(result["solution"]["storage"]["1"]["ps"], -0.0596443; atol=1e-3))
        @test all(isapprox.(result["solution"]["storage"]["2"]["ps"], -0.0793700; atol=1e-3))
    end

    @testset "5-bus storage nfa opf_strg" begin
        result = run_mc_opf(mp_data, NFAPowerModel, ipopt_solver; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 43169.9; atol = 1e0)

        @test isapprox(result["solution"]["storage"]["1"]["se"],  0.0; atol = 1e0)
        @test isapprox(result["solution"]["storage"]["2"]["se"],  0.0; atol = 1e0)

        @test all(isapprox.(result["solution"]["storage"]["1"]["ps"], -0.0596443; atol=1e-3))
        @test all(isapprox.(result["solution"]["storage"]["2"]["ps"], -0.0793700; atol=1e-3))
    end
end
