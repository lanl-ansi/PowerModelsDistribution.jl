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

        @test all(isapprox.(sol["solution"]["bus"]["2"]["vm"].values, 0.984406; atol=1e-4))
        @test all(isapprox.(sol["solution"]["bus"]["2"]["va"].values, [-0.0114755, -2.10587, 2.08292]; atol=deg2rad(0.2)))

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
    end

    @testset "3-bus balanced" begin
        tppm = TPPMs.parse_file("../test/data/opendss/case3_balanced.dss")
        sol = TPPMs.run_tp_pf(tppm, PMs.ACPPowerModel, ipopt_solver)

        @test sol["status"] == :LocalOptimal

        for (bus, va, vm) in zip(["1", "2", "3"], [[0.0, -2.0944, 2.0944], [-0.00057819, -2.09497, 2.09382], [-0.0011636, -2.09556, 2.09323]], [0.9959, 0.986976, 0.976611])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"].values, va; atol=deg2rad(0.2)))
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
                                [[0.0, -2.0944, 2.0944], [-0.00391758, -2.09638, 2.09654], [-0.00845653, -2.09863, 2.09917]],
                                [0.9959, [0.980939, 0.989362, 0.987041], [0.963549, 0.98176, 0.976782]])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"].values, va; atol=deg2rad(0.2)))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"].values, vm; atol=2e-3))
        end

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214835; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00932693; atol=1e-4)
    end
end
