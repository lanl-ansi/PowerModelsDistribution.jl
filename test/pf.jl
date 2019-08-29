@testset "test ac polar pf" begin
    @testset "5-bus independent meshed network" begin
        result = run_ac_mc_pf("../test/data/matlab/case5_i_m_b.m", ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "5-bus coupled meshed network (a)" begin
        result = run_ac_mc_pf("../test/data/matlab/case5_c_m_a.m", ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "5-bus coupled meshed network (b)" begin
        result = run_ac_mc_pf("../test/data/matlab/case5_c_m_b.m", ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "5-bus independent radial w/ shunts" begin
        result = run_ac_mc_pf("../test/data/matlab/case5_i_r_b.m", ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "2-bus diagonal" begin
        pmd = PMD.parse_file("../test/data/opendss/case2_diag.dss")
        sol = PMD.run_mc_pf(pmd, PMs.ACPPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        @test all(isapprox.(sol["solution"]["bus"]["2"]["vm"].values, 0.984377; atol=1e-4))
        @test all(isapprox.(sol["solution"]["bus"]["2"]["va"].values, [PMD._wrap_to_pi(2*pi/pmd["conductors"]*(1-c) - deg2rad(0.79)) for c in 1:pmd["conductors"]]; atol=deg2rad(0.2)))

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
    end

    @testset "3-bus balanced" begin
        pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
        sol = PMD.run_mc_pf(pmd, PMs.ACPPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        for (bus, va, vm) in zip(["1", "2", "3"], [0.0, deg2rad(-0.08), deg2rad(-0.17)], [0.9959, 0.986976, 0.976611])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"].values, [PMD._wrap_to_pi(2*pi/pmd["conductors"]*(1-c) .+ va) for c in 1:pmd["conductors"]]; atol=deg2rad(0.2)))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"].values, vm; atol=1e-3))
        end

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183456; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923328; atol=1e-4)

        @testset "3-bus balanced no linecode basefreq defined" begin
            pmd2 = PMD.parse_file("../test/data/opendss/case3_balanced_basefreq.dss")
            sol2 = PMD.run_mc_pf(pmd2, PMs.ACPPowerModel, ipopt_solver)

            @test all(all(isapprox.(bus["vm"].values, sol2["solution"]["bus"][i]["vm"]; atol=1e-8)) for (i, bus) in sol["solution"]["bus"])
            @test all(all(isapprox.(bus["va"].values, sol2["solution"]["bus"][i]["va"]; atol=1e-8)) for (i, bus) in sol["solution"]["bus"])
            @test all(isapprox(sum(sol["solution"]["gen"]["1"][field] * sol["solution"]["baseMVA"]), sum(sol2["solution"]["gen"]["1"][field] * sol2["solution"]["baseMVA"]); atol=1e-8) for field in ["pg", "qg"])
        end
    end

    @testset "3-bus unbalanced" begin
        pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced.dss")
        sol = PMD.run_mc_pf(pmd, PMs.ACPPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        for (bus, va, vm) in zip(["1", "2", "3"],
                                [0.0, deg2rad.([-0.22, -0.11, 0.12]), deg2rad.([-0.48, -0.24, 0.27])],
                                [0.9959, [0.98094, 0.989365, 0.987043], [0.96355, 0.981767, 0.976786]])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"].values, PMD._wrap_to_pi([2*pi/pmd["conductors"]*(1-c) for c in 1:pmd["conductors"]] .+ va); atol=deg2rad(0.2)))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"].values, vm; atol=2e-3))
        end

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214835; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00932693; atol=1e-4)
    end

    @testset "5-bus 3-phase ac pf case" begin
        mp_data = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss")
        result = run_mc_pf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol = 1e-4)

        @test isapprox(result["solution"]["gen"]["1"]["pg"][1], 0.00015328882711864364; atol = 1e-5)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][2], 0.0001993266190368713; atol = 1e-5)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][3], 0.0002480554356591965; atol = 1e-5)

        @test isapprox(result["solution"]["bus"]["2"]["vm"][1], 0.9733455037213229; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][2], 0.9647241898338335; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][3], 0.9555739064829893; atol = 1e-3)
    end

    @testset "matrix branch shunts" begin
        sol = PMD.run_ac_mc_pf_lm("../test/data/opendss/case_mxshunt.dss", ipopt_solver)

        # these results were obtained from OpenDSS; largest mismatch was 5E-9
        isapprox(sol["solution"]["bus"]["2"]["vm"][1], 0.9873988561202298, atol=1E-8)
        isapprox(sol["solution"]["bus"]["2"]["vm"][2], 0.9813000619074207, atol=1E-8)
        isapprox(sol["solution"]["bus"]["2"]["vm"][3], 1.0035368353626686, atol=1E-8)
    end
end

@testset "test ac rectangular pf" begin
    @testset "5-bus independent meshed network" begin
        result = PMD.run_mc_pf("../test/data/matlab/case5_i_m_b.m", PMs.ACRPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "5-bus coupled meshed network (a)" begin
        result = PMD.run_mc_pf("../test/data/matlab/case5_c_m_a.m", PMs.ACRPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "5-bus coupled meshed network (b)" begin
        result = PMD.run_mc_pf("../test/data/matlab/case5_c_m_b.m", PMs.ACRPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "5-bus independent radial w/ shunts" begin
        result = PMD.run_mc_pf("../test/data/matlab/case5_i_r_b.m", PMs.ACRPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end

    @testset "2-bus diagonal" begin
        pmd = PMD.parse_file("../test/data/opendss/case2_diag.dss")
        sol = PMD.run_mc_pf(pmd, PMs.ACRPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        @test all(isapprox.(sol["solution"]["bus"]["2"]["vm"].values, 0.984377; atol=1e-4))
        @test all(isapprox.(sol["solution"]["bus"]["2"]["va"].values, [PMD._wrap_to_pi(2*pi/pmd["conductors"]*(1-c) - deg2rad(0.79)) for c in 1:pmd["conductors"]]; atol=deg2rad(0.2)))

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
    end

    @testset "3-bus balanced" begin
        pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
        sol = PMD.run_mc_pf(pmd, PMs.ACRPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        for (bus, va, vm) in zip(["1", "2", "3"], [0.0, deg2rad(-0.08), deg2rad(-0.17)], [0.9959, 0.986976, 0.976611])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"].values, [PMD._wrap_to_pi(2*pi/pmd["conductors"]*(1-c) .+ va) for c in 1:pmd["conductors"]]; atol=deg2rad(0.2)))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"].values, vm; atol=1e-3))
        end

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183456; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923328; atol=1e-4)
    end

    @testset "3-bus unbalanced" begin
        pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced.dss")
        sol = PMD.run_mc_pf(pmd, PMs.ACRPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        for (bus, va, vm) in zip(["1", "2", "3"],
                                [0.0, deg2rad.([-0.22, -0.11, 0.12]), deg2rad.([-0.48, -0.24, 0.27])],
                                [0.9959, [0.98094, 0.989365, 0.987043], [0.96355, 0.981767, 0.976786]])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"].values, PMD._wrap_to_pi([2*pi/pmd["conductors"]*(1-c) for c in 1:pmd["conductors"]] .+ va); atol=deg2rad(0.2)))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"].values, vm; atol=2e-3))
        end

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214835; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00932693; atol=1e-4)
    end

    @testset "5-bus 3-phase ac pf case" begin
        mp_data = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss")
        result = run_mc_pf(mp_data, PMs.ACRPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol = 1e-4)

        @test isapprox(result["solution"]["gen"]["1"]["pg"][1], 0.00015328882711864364; atol = 1e-5)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][2], 0.0001993266190368713; atol = 1e-5)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][3], 0.0002480554356591965; atol = 1e-5)

        @test isapprox(result["solution"]["bus"]["2"]["vm"][1], 0.9733455037213229; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][2], 0.9647241898338335; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][3], 0.9555739064829893; atol = 1e-3)
    end

    @testset "matrix branch shunts" begin
        data_pmd = PMD.parse_file("../test/data/opendss/case_mxshunt.dss")
        pm = PMs.build_model(data_pmd, PMs.ACRPowerModel, PMD.post_mc_pf_lm, ref_extensions=[PMD.ref_add_arcs_trans!], multiconductor=true)
        sol = PMs.optimize_model!(pm, ipopt_solver)

        # these results were obtained from OpenDSS; largest mismatch was 5E-9
        isapprox(sol["solution"]["bus"]["2"]["vm"][1], 0.9873988561202298, atol=1E-8)
        isapprox(sol["solution"]["bus"]["2"]["vm"][2], 0.9813000619074207, atol=1E-8)
        isapprox(sol["solution"]["bus"]["2"]["vm"][3], 1.0035368353626686, atol=1E-8)
    end
end

@testset "test soc pf" begin
    @testset "5-bus coupled meshed network (b)" begin
        result = run_mc_pf("../test/data/matlab/case5_c_m_b.m", PMs.SOCWRPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol=1e-2)
    end
end
