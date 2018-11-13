@testset "test ac polar pf" begin
    @testset "5-bus independent meshed network" begin
        result = run_ac_tp_pf("../test/data/matlab/case5_i_m_b.m", ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "5-bus coupled meshed network (a)" begin
        result = run_ac_tp_pf("../test/data/matlab/case5_c_m_a.m", ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "5-bus coupled meshed network (b)" begin
        result = run_ac_tp_pf("../test/data/matlab/case5_c_m_b.m", ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "5-bus independent radial w/ shunts" begin
        result = run_ac_tp_pf("../test/data/matlab/case5_i_r_b.m", ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "2-bus diagonal" begin
        tppm = TPPMs.parse_file("../test/data/opendss/case2_diag.dss")
        sol = TPPMs.run_tp_pf(tppm, PMs.ACPPowerModel, ipopt_solver)

        @test sol["status"] == :LocalOptimal

        @test all(isapprox.(sol["solution"]["bus"]["2"]["vm"].values, 0.984377; atol=1e-4))
        @test all(isapprox.(sol["solution"]["bus"]["2"]["va"].values, [TPPMs.wraptopi(2*pi/tppm["conductors"]*(1-c) - deg2rad(0.79)) for c in 1:tppm["conductors"]]; atol=deg2rad(0.2)))

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
    end

    @testset "3-bus balanced" begin
        tppm = TPPMs.parse_file("../test/data/opendss/case3_balanced.dss")
        sol = TPPMs.run_tp_pf(tppm, PMs.ACPPowerModel, ipopt_solver)

        @test sol["status"] == :LocalOptimal

        for (bus, va, vm) in zip(["1", "2", "3"], [0.0, deg2rad(-0.08), deg2rad(-0.17)], [0.9959, 0.986976, 0.976611])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"].values, [TPPMs.wraptopi(2*pi/tppm["conductors"]*(1-c) .+ va) for c in 1:tppm["conductors"]]; atol=deg2rad(0.2)))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"].values, vm; atol=1e-3))
        end

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183456; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923328; atol=1e-4)
    end

    @testset "3-bus unbalanced" begin
        tppm = TPPMs.parse_file("../test/data/opendss/case3_unbalanced.dss")
        sol = TPPMs.run_tp_pf(tppm, PMs.ACPPowerModel, ipopt_solver)

        @test sol["status"] == :LocalOptimal

        for (bus, va, vm) in zip(["1", "2", "3"],
                                [0.0, deg2rad.([-0.22, -0.11, 0.12]), deg2rad.([-0.48, -0.24, 0.27])],
                                [0.9959, [0.98094, 0.989365, 0.987043], [0.96355, 0.981767, 0.976786]])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"].values, TPPMs.wraptopi([2*pi/tppm["conductors"]*(1-c) for c in 1:tppm["conductors"]] .+ va); atol=deg2rad(0.2)))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"].values, vm; atol=2e-3))
        end

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214835; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00932693; atol=1e-4)
    end

    @testset "5-bus 3-phase ac pf case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/case5_phase_drop.dss")
        result = run_tp_pf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol = 1e-4)

        @test isapprox(result["solution"]["gen"]["1"]["pg"][1], 0.00015328882711864364; atol = 1e-5)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][2], 0.0001993266190368713; atol = 1e-5)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][3], 0.0002480554356591965; atol = 1e-5)

        @test isapprox(result["solution"]["bus"]["2"]["vm"][1], 0.9733455037213229; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][2], 0.9647241898338335; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][3], 0.9555739064829893; atol = 1e-3)
    end
end


@testset "test soc pf" begin
    @testset "5-bus coupled meshed network (b)" begin
        result = run_tp_pf("../test/data/matlab/case5_c_m_b.m", PMs.SOCWRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end
end
