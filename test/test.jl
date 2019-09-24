@info "running storage optimal power flow tests"

@testset "test storage opf" begin
    @testset "5-bus storage acp opf_strg" begin
        mp_data = PowerModels.parse_file("../test/data/matpower/case5_strg.m")
        PMD.make_multiconductor!(mp_data, 3)

        result = PMD.run_mc_opf(mp_data, PowerModels.ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 52299.2; atol = 1e0)

        @test isapprox(result["solution"]["storage"]["1"]["se"],  0.0; atol = 1e0)
        @test isapprox(result["solution"]["storage"]["2"]["se"],  0.0; atol = 1e0)
        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["storage"]["1"]["ps"][c], -0.0596928; atol = 1e-3)
            @test isapprox(result["solution"]["storage"]["2"]["ps"][c], -0.0794600; atol = 1e-3)
        end
    end

    @testset "5-bus storage dcp opf_strg" begin
        mp_data = PowerModels.parse_file("../test/data/matpower/case5_strg.m")
        PMD.make_multiconductor!(mp_data, 3)

        result = PMD.run_mc_opf(mp_data, PowerModels.DCPPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 52059.6; atol = 1e0)

        @test isapprox(result["solution"]["storage"]["1"]["se"],  0.0; atol = 1e0)
        @test isapprox(result["solution"]["storage"]["2"]["se"],  0.0; atol = 1e0)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["storage"]["1"]["ps"][c], -0.0596443; atol = 1e-3)
            @test isapprox(result["solution"]["storage"]["2"]["ps"][c], -0.0793700; atol = 1e-3)
        end
    end

    @testset "5-bus storage nfa opf_strg" begin
        mp_data = PowerModels.parse_file("../test/data/matpower/case5_strg.m")
        PMD.make_multiconductor!(mp_data, 3)

        result = PMD.run_mc_opf(mp_data, PowerModels.NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 43169.9; atol = 1e0)

        @test isapprox(result["solution"]["storage"]["1"]["se"],  0.0; atol = 1e0)
        @test isapprox(result["solution"]["storage"]["2"]["se"],  0.0; atol = 1e0)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["storage"]["1"]["ps"][c], -0.0596443; atol = 1e-3)
            @test isapprox(result["solution"]["storage"]["2"]["ps"][c], -0.0793700; atol = 1e-3)
        end
    end
end
