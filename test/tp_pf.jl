@testset "test ac polar pf" begin
    @testset "5-bus independent meshed network" begin
        result = run_ac_tp_pf("../test/data/matlab/case5_i_m_b.m", ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "5-bus coupled meshed network" begin
        result = run_ac_tp_pf("../test/data/matlab/case5_c_m_a.m", ipopt_solver)

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
        @test all(isapprox.(sol["solution"]["bus"]["2"]["va"].values, [2 * pi / tppm["conductors"] * (c - 1) - deg2rad(0.79) for c in 1:tppm["conductors"]]; atol=deg2rad(0.2)))

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
    end

    @testset "3-bus balanced" begin
        tppm = TPPMs.parse_file("../test/data/opendss/case3_balanced.dss")
        sol = TPPMs.run_tp_pf(tppm, PMs.ACPPowerModel, ipopt_solver)

        @test sol["status"] == :LocalOptimal

        for (bus, va, vm) in zip(["1", "2", "3"], [0.0, deg2rad(-0.08), deg2rad(-0.17)], [0.9959, 0.986559, 0.97572])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"].values, [2 * pi / tppm["conductors"] * (c - 1) + va for c in 1:tppm["conductors"]]; atol=deg2rad(0.2)))
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
                                [0.0, deg2rad.([-0.30, 0.09, -0.17]), deg2rad.([-0.65, 0.20, -0.36])],
                                [0.9959, [0.980269, 0.986645, 0.989161], [0.962159, 0.975897, 0.981341]])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"].values, [2 * pi / tppm["conductors"] * (c - 1) for c in 1:tppm["conductors"]] + va; atol=deg2rad(0.2)))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"].values, vm; atol=2e-3))
        end

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214835; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00932693; atol=1e-4)
    end
end