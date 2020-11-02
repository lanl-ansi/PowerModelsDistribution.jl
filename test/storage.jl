@info "running storage tests"

@testset "test storage opf" begin
    mp_data = PM.parse_file("../test/data/matpower/case5_strg.m")
    make_multiconductor!(mp_data, 3)

    @testset "5-bus storage acp opf" begin
        result = run_mc_opf(mp_data, ACPPowerModel, ipopt_solver; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 52299.2; atol = 1e0)

        @test isapprox(result["solution"]["storage"]["1"]["se"],  0.0; atol=1e0)
        @test isapprox(result["solution"]["storage"]["2"]["se"],  0.0; atol=1e0)

        @test all(isapprox.(result["solution"]["storage"]["1"]["ps"], -0.0596928; atol=1e-3))
        @test all(isapprox.(result["solution"]["storage"]["2"]["ps"], -0.0794600; atol=1e-3))
    end

    @testset "5-bus storage arc opf" begin
        result = run_mc_opf(mp_data, ACRPowerModel, ipopt_solver; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 52299.2; atol = 1e0)

        @test isapprox(result["solution"]["storage"]["1"]["se"],  0.0; atol=1e0)
        @test isapprox(result["solution"]["storage"]["2"]["se"],  0.0; atol=1e0)

        @test all(isapprox.(result["solution"]["storage"]["1"]["ps"], -0.0596928; atol=1e-3))
        @test all(isapprox.(result["solution"]["storage"]["2"]["ps"], -0.0794600; atol=1e-3))
    end

    @testset "5-bus storage dcp opf" begin
        result = run_mc_opf(mp_data, DCPPowerModel, ipopt_solver; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 52059.6; atol = 1e0)

        @test isapprox(result["solution"]["storage"]["1"]["se"],  0.0; atol = 1e0)
        @test isapprox(result["solution"]["storage"]["2"]["se"],  0.0; atol = 1e0)

        @test all(isapprox.(result["solution"]["storage"]["1"]["ps"], -0.0596443; atol=1e-3))
        @test all(isapprox.(result["solution"]["storage"]["2"]["ps"], -0.0793700; atol=1e-3))
    end

    @testset "5-bus storage lpubfdiag opf" begin
        result = run_mc_opf(mp_data, LPUBFDiagPowerModel, ipopt_solver; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 43177.3; atol = 1e0)

        @test isapprox(result["solution"]["storage"]["1"]["se"],  0.0; atol = 1e0)
        @test isapprox(result["solution"]["storage"]["2"]["se"],  0.0; atol = 1e0)

        @test all(isapprox.(result["solution"]["storage"]["1"]["ps"], -0.059; atol=1e-3))
        @test all(isapprox.(result["solution"]["storage"]["2"]["ps"], -0.079; atol=1e-3))
    end

    @testset "5-bus storage nfa opf" begin
        result = run_mc_opf(mp_data, NFAPowerModel, ipopt_solver; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 43169.9; atol = 1e0)

        @test isapprox(result["solution"]["storage"]["1"]["se"],  0.0; atol = 1e0)
        @test isapprox(result["solution"]["storage"]["2"]["se"],  0.0; atol = 1e0)

        @test all(isapprox.(result["solution"]["storage"]["1"]["ps"], -0.0596443; atol=1e-3))
        @test all(isapprox.(result["solution"]["storage"]["2"]["ps"], -0.0793700; atol=1e-3))
    end
end


@testset "test storage pf" begin
    eng = parse_file("../test/data/opendss/case3_balanced_battery.dss")

    @testset "3-bus balanced battery acp pf" begin
        result = run_mc_pf(eng, ACPPowerModel, ipopt_solver; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"], 0.98697; atol=1e-5))
    end
    @testset "3-bus balanced battery acr pf" begin
        result = run_mc_pf(eng, ACRPowerModel, ipopt_solver; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test all(isapprox.(calc_vm_acr(result, "primary"), 0.98697; atol=1e-5))
    end
    # TODO this test is unstable on travis, fails randomly, needs better test
    # @testset "3-bus balanced battery lpubfdiag pf" begin
    #     result = run_mc_pf(eng, LPUBFDiagPowerModel, ipopt_solver; make_si=false)

    #     @test result["termination_status"] == LOCALLY_SOLVED
    #     @test all(isapprox.(result["solution"]["bus"]["primary"]["w"], 0.99767; atol=1e-2))
    # end
end

@testset "test storage mld" begin
    mp_data = PM.parse_file("../test/data/matpower/case5_mld_strg.m")
    make_multiconductor!(mp_data, 3)

    @testset "5-bus mld storage acp mld" begin
        result = run_mc_mld(mp_data, ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 137.17; atol = 1e-2)

        @test all(isapprox(gen["gen_status"], 1.0; atol=1e-6) for (_,gen) in result["solution"]["gen"])
        @test isapprox(sum(sum(load["pd"]) for (_,load) in result["solution"]["load"]), 18.09; atol=1e-2)
    end

    @testset "5-bus mld storage acr mld" begin
        result = run_mc_mld(mp_data, ACRPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 137.17; atol = 1e-2)

        @test all(isapprox(gen["gen_status"], 1.0; atol=1e-6) for (_,gen) in result["solution"]["gen"])
        @test isapprox(sum(sum(load["pd"]) for (_,load) in result["solution"]["load"]), 18.09; atol=1e-2)
    end

    @testset "5-bus mld storage lpubfdiag mld" begin
        result = run_mc_mld(mp_data, LPUBFDiagPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 136.94; atol = 1e-1)

        @test all(isapprox(gen["gen_status"], 1.0; atol=1e-6) for (_,gen) in result["solution"]["gen"])
        @test isapprox(sum(sum(load["pd"]) for (_,load) in result["solution"]["load"]), 18.11; atol=1e-2)
    end

    @testset "5-bus mld storage nfa mld" begin
        result = run_mc_mld(mp_data, NFAPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 131.7; atol = 1e-1)

        @test all(isapprox(gen["gen_status"], 1.0; atol=1e-6) for (_,gen) in result["solution"]["gen"])
        @test isapprox(sum(sum(load["pd"]) for (_,load) in result["solution"]["load"]), 18.69; atol=1e-2)
    end
end
