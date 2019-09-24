@info "running power flow (pf) tests"

@testset "test pf" begin
    @testset "test pmd matlab pf" begin
        @testset "5-bus independent meshed network acp pf" begin
            result = run_ac_mc_pf("../test/data/matlab/case5_i_m_b.m", ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0; atol=1e-2)
        end

        @testset "5-bus independent meshed network acr pf" begin
            result = PMD.run_mc_pf("../test/data/matlab/case5_i_m_b.m", PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0; atol=1e-2)
        end

        @testset "5-bus coupled meshed network (a) acp pf" begin
            result = run_ac_mc_pf("../test/data/matlab/case5_c_m_a.m", ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0; atol=1e-2)
        end

        @testset "5-bus coupled meshed network (a) acr pf" begin
            result = PMD.run_mc_pf("../test/data/matlab/case5_c_m_a.m", PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0; atol=1e-2)
        end

        @testset "5-bus coupled meshed network (b) acp pf" begin
            result = run_ac_mc_pf("../test/data/matlab/case5_c_m_b.m", ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0; atol=1e-2)
        end

        @testset "5-bus coupled meshed network (b) acr pf" begin
            result = PMD.run_mc_pf("../test/data/matlab/case5_c_m_b.m", PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0; atol=1e-2)
        end

        @testset "5-bus coupled meshed network (b) soc pf" begin
            result = run_mc_pf("../test/data/matlab/case5_c_m_b.m", PMs.SOCWRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0; atol=1e-2)
        end

        @testset "5-bus independent radial w/ shunts acp pf" begin
            result = run_ac_mc_pf("../test/data/matlab/case5_i_r_b.m", ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0; atol=1e-2)
        end

        @testset "5-bus independent radial w/ shunts acr pf" begin
            result = PMD.run_mc_pf("../test/data/matlab/case5_i_r_b.m", PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0; atol=1e-2)
        end
    end

    @testset "test opendss pf" begin
        @testset "2-bus diagonal acp pf" begin
            pmd = PMD.parse_file("../test/data/opendss/case2_diag.dss")
            sol = PMD.run_mc_pf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            @test all(isapprox.(sol["solution"]["bus"]["2"]["vm"].values, 0.984377; atol=1e-4))
            @test all(isapprox.(sol["solution"]["bus"]["2"]["va"].values, [PMD._wrap_to_pi(2*pi/pmd["conductors"]*(1-c) - deg2rad(0.79)) for c in 1:pmd["conductors"]]; atol=deg2rad(0.2)))

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
        end

        @testset "2-bus diagonal acr pf" begin
            pmd = PMD.parse_file("../test/data/opendss/case2_diag.dss")
            sol = PMD.run_mc_pf(pmd, PMs.ACRPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            @test all(isapprox.(sol["solution"]["bus"]["2"]["vm"].values, 0.984377; atol=1e-4))
            @test all(isapprox.(sol["solution"]["bus"]["2"]["va"].values, [PMD._wrap_to_pi(2*pi/pmd["conductors"]*(1-c) - deg2rad(0.79)) for c in 1:pmd["conductors"]]; atol=deg2rad(0.2)))

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
        end

        @testset "3-bus balanced acp pf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
            sol = PMD.run_mc_pf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            for (bus, va, vm) in zip(["1", "2", "3"], [0.0, deg2rad(-0.08), deg2rad(-0.17)], [0.9959, 0.986976, 0.976611])
                @test all(isapprox.(sol["solution"]["bus"][bus]["va"].values, [PMD._wrap_to_pi(2*pi/pmd["conductors"]*(1-c) .+ va) for c in 1:pmd["conductors"]]; atol=deg2rad(0.2)))
                @test all(isapprox.(sol["solution"]["bus"][bus]["vm"].values, vm; atol=1e-3))
            end

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183456; atol=1e-5)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923328; atol=1e-4)
        end

        @testset "3-bus balanced acr pf" begin
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

        @testset "3-bus balanced no linecode basefreq defined acp pf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
            sol = PMD.run_mc_pf(pmd, PMs.ACPPowerModel, ipopt_solver)

            pmd2 = PMD.parse_file("../test/data/opendss/case3_balanced_basefreq.dss")
            sol2 = PMD.run_mc_pf(pmd2, PMs.ACPPowerModel, ipopt_solver)

            @test all(all(isapprox.(bus["vm"].values, sol2["solution"]["bus"][i]["vm"]; atol=1e-8)) for (i, bus) in sol["solution"]["bus"])
            @test all(all(isapprox.(bus["va"].values, sol2["solution"]["bus"][i]["va"]; atol=1e-8)) for (i, bus) in sol["solution"]["bus"])
            @test all(isapprox(sum(sol["solution"]["gen"]["1"][field] * sol["solution"]["baseMVA"]), sum(sol2["solution"]["gen"]["1"][field] * sol2["solution"]["baseMVA"]); atol=1e-8) for field in ["pg", "qg"])
        end

        @testset "3-bus balanced w/ switch acp pf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced_switch.dss")
            sol = PMD.run_mc_pf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            for (bus, va, vm) in zip(["1", "2", "3"], [0.0, 0.0, deg2rad(-0.04)], [0.9959, 0.995729, 0.985454])
                @test all(isapprox.(sol["solution"]["bus"][bus]["va"].values, [PMD._wrap_to_pi(2*pi/pmd["conductors"]*(1-c) .+ va) for c in 1:pmd["conductors"]]; atol=deg2rad(0.2)))
                @test all(isapprox.(sol["solution"]["bus"][bus]["vm"].values, vm; atol=1e-3))
            end
        end

        @testset "3-bus unbalanced acp pf" begin
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

        @testset "3-bus unbalanced acr pf" begin
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

        @testset "3-bus unbalanced w/ assymetric linecode & phase order swap acp pf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced_assym_swap.dss")
            sol = PMD.run_ac_mc_pf(pmd, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            @test all(isapprox.(sol["solution"]["bus"]["2"]["vm"], [0.983453, 0.98718, 0.981602]; atol=1e-5))
            @test all(isapprox.(sol["solution"]["bus"]["2"]["va"], deg2rad.([-0.07, -120.19, 120.29]); atol=1e-2))
        end

        @testset "5-bus phase drop acp pf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss")
            result = run_mc_pf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0; atol = 1e-4)

            @test all(isapprox.(result["solution"]["bus"]["2"]["vm"], [0.973519, 0.964902, 0.956465]; atol = 1e-3))
        end

        @testset "5-bus phase drop acr pf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss")
            result = run_mc_pf(mp_data, PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0; atol = 1e-4)

            @test isapprox(result["solution"]["bus"]["2"]["vm"][1], 0.973519; atol = 1e-4)
            @test isapprox(result["solution"]["bus"]["2"]["vm"][2], 0.964902; atol = 1e-4)
            @test isapprox(result["solution"]["bus"]["2"]["vm"][3], 0.956465; atol = 1e-4)
        end

        @testset "matrix branch shunts acp pf" begin
            sol = PMD.run_ac_mc_pf_lm("../test/data/opendss/case_mxshunt.dss", ipopt_solver)

            @test all(isapprox.(sol["solution"]["bus"]["2"]["vm"], [0.987399, 0.981300, 1.003536]; atol=1E-6))
        end

        @testset "matrix branch shunts acr pf" begin
            data_pmd = PMD.parse_file("../test/data/opendss/case_mxshunt.dss")
            pm = PMs.build_model(data_pmd, PMs.ACRPowerModel, PMD.post_mc_pf_lm, ref_extensions=[PMD.ref_add_arcs_trans!], multiconductor=true)
            sol = PMs.optimize_model!(pm, ipopt_solver)

            @test all(isapprox.(sol["solution"]["bus"]["2"]["vm"], [0.987399, 0.981299, 1.003537]; atol=1E-6))
        end

        @testset "virtual sourcebus creation acp pf" begin
            result = run_ac_mc_pf("../test/data/opendss/virtual_sourcebus.dss", ipopt_solver)
            @test result["termination_status"] == PMs.LOCALLY_SOLVED

            @test all(all(isapprox.(result["solution"]["bus"]["$n"]["vm"].values, [0.961352, 0.999418, 1.00113]; atol=1e-6)) for n in [1, 2])
            @test all(all(isapprox.(rad2deg.(result["solution"]["bus"]["$n"]["va"].values), [-1.25, -120.06, 120.0]; atol=1e-1)) for n in [1, 2])
        end
    end
end

@testset "test pf_bf" begin
    @testset "test opendss pf_bf" begin
        @testset "3-bus unbalanced acp pf_bf branch_flows original_variables" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced.dss")
            sol = PMD.run_mc_pf_bf(pmd, LPLinUBFPowerModel, ipopt_solver, setting=Dict("output"=>Dict("branch_flows"=>true, "original_variables"=>true)))

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.0209; atol=1e-4)
            @test all(isapprox.(sol["solution"]["bus"]["1"]["vm"], 0.9959; atol=1e-4))
            @test all(haskey(sol["solution"]["branch"]["1"],key) for key in [ "qf_ut" "pt" "qt_lt" "cci" "ccr" "ctm" "cta" "qt_ut" "qf_lt" "pt_lt" "pf_lt" "ccm" "cfm" "cc" "qf" "pt_ut" "qt" "cfa" "pf" "pf_ut" ])
        end
    end
end
