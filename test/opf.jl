@info "running optimal power flow (opf) tests"

@testset "test opf" begin
    @testset "test matpower opf" begin
        @testset "5-bus matpower acp opf" begin
            mp_data = PMs.parse_file("../test/data/matpower/case5.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 45522.096; atol=1e-1)

            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.3999999; atol=1e-3) for c in 1:mp_data["conductors"])
            @test all(isapprox(result["solution"]["bus"]["2"]["va"][c], PMD._wrap_to_pi(-0.0538204+2*pi/mp_data["conductors"]*(1-c)); atol=1e-5) for c in 1:mp_data["conductors"])
        end

        @testset "5-bus matpower acr opf" begin
            mp_data = PMs.parse_file("../test/data/matpower/case5.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 45522.096; atol=1e-1)

            calc_va(id) = atan.(result["solution"]["bus"][id]["vi"], result["solution"]["bus"][id]["vr"])
            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.3999999; atol=1e-3) for c in 1:mp_data["conductors"])
            @test all(isapprox(calc_va("2")[c], PMD._wrap_to_pi(-0.0538204+2*pi/mp_data["conductors"]*(1-c)); atol=1e-5) for c in 1:mp_data["conductors"])
        end

        # @testset "5-bus matpower soc opf" begin
        #     mp_data = PMs.parse_file("../test/data/matpower/case5.m")
        #     PMD.make_multiconductor!(mp_data, 3)
        #     result = run_mc_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)
        #
        #     @test result["termination_status"] == PMs.LOCALLY_SOLVED
        #     @test isapprox(result["objective"], 45365.17; atol=1e-1)
        #
        #     @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.3999999; atol=1e-3) for c in 1:mp_data["conductors"])
        # end

        @testset "30-bus matpower acp opf" begin
            mp_data = PMs.parse_file("../test/data/matpower/case30.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 614.007; atol=1e-1)

            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  2.192189; atol=1e-3) for c in 1:mp_data["conductors"])
            @test all(isapprox(result["solution"]["bus"]["2"]["va"][c], PMD._wrap_to_pi(-0.071853+2*pi/mp_data["conductors"]*(1-c)); atol=1e-4) for c in 1:mp_data["conductors"])
        end

        @testset "30-bus matpower acr opf" begin
            mp_data = PMs.parse_file("../test/data/matpower/case30.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 614.007; atol=1e-1)

            calc_va(id) = atan.(result["solution"]["bus"][id]["vi"], result["solution"]["bus"][id]["vr"])
            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  2.192189; atol=1e-3) for c in 1:mp_data["conductors"])
            @test all(isapprox(calc_va("2")[c], PMD._wrap_to_pi(-0.071853+2*pi/mp_data["conductors"]*(1-c)); atol=1e-4) for c in 1:mp_data["conductors"])
        end

        @testset "30-bus matpower dcp opf" begin
            mp_data = PMs.parse_file("../test/data/matpower/case30.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf(mp_data, PMs.DCPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 566.112; atol=1e-1)
        end

        @testset "30-bus matpower nfa opf" begin
            mp_data = PMs.parse_file("../test/data/matpower/case30.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf(mp_data, PMs.NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 458.006; atol=1e-1)
        end

        # @testset "30-bus matpower soc opf" begin
        #     mp_data = PMs.parse_file("../test/data/matpower/case30.m")
        #     PMD.make_multiconductor!(mp_data, 3)
        #     result = run_mc_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)
        #
        #     @test result["termination_status"] == PMs.LOCALLY_SOLVED
        #     @test isapprox(result["objective"], 517.588; atol=1e-1)
        #
        #     @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  2.821313; atol=1e-3) for c in 1:mp_data["conductors"])
        # end
    end

    @testset "test pmd matlab opf" begin
        @testset "5-bus independent radial identical acp opf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_i_r_a.m")
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 55451.7; atol=1e-1)

            @test all(isapprox(result["solution"]["gen"]["1"]["qg"][c], 0.039742; atol=1e-4) for c in 1:mp_data["conductors"])
            @test all(isapprox(result["solution"]["bus"]["2"]["va"][c], PMD._wrap_to_pi(0.048896+2*pi/mp_data["conductors"]*(1-c)); atol=1e-4) for c in 1:mp_data["conductors"])
        end

        @testset "5-bus independent radial different acp opf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_i_r_b.m")
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 56091.3; atol=1e-1)

            @test isapprox(result["solution"]["gen"]["1"]["qg"][1],  0.105276; atol=1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][1],  0.0575114; atol=1e-3)

            @test all(isapprox(result["solution"]["gen"]["1"]["qg"][c],  0.0897773; atol=1e-3) for c in 2:mp_data["conductors"])
            @test all(isapprox(result["solution"]["bus"]["2"]["va"][c],  PMD._wrap_to_pi(0.052544+2*pi/mp_data["conductors"]*(1-c)); atol=1e-3) for c in 2:mp_data["conductors"])
        end

        @testset "5-bus independent meshed different acp opf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_i_m_b.m")
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 52964.4; atol=1e-1)

            @test all(isapprox(result["solution"]["gen"]["1"]["qg"][c],  0.3; atol=1e-3) for c in 1:mp_data["conductors"])
            @test all(isapprox(result["solution"]["bus"]["2"]["va"][c], PMD._wrap_to_pi(-0.0135651+2*pi/mp_data["conductors"]*(1-c)); atol=1e-3) for c in 1:mp_data["conductors"])
        end

        @testset "5-bus coupled meshed acp opf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 53272.9; atol=1e-1)

            @test all(isapprox.(result["solution"]["gen"]["1"]["qg"], 0.3; atol=1e-3))

            @test all(isapprox.(result["solution"]["bus"]["2"]["va"], [-0.0139580, -2.1069476, 2.0808321]; atol=1e-3))
        end

        @testset "5-bus coupled meshed acr opf" begin
            mp_data = PowerModelsDistribution.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_mc_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 53272.9; atol=1e-1)
        end

        @testset "5-bus coupled meshed dcp opf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_mc_opf(mp_data, PMs.DCPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 55640.2; atol=1e-1)
        end

        @testset "5-bus coupled meshed nfa opf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_mc_opf(mp_data, PMs.NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 44700.0; atol=1e-1)
        end

        # @testset "5-bus coupled meshed soc opf" begin
        #     mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_a.m")
        #     result = run_mc_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)
        #
        #     @test result["termination_status"] == PMs.LOCALLY_SOLVED
        #     @test isapprox(result["objective"], -0.000272; atol=1e-3)
        #
        #     @test all(isapprox.(result["solution"]["gen"]["1"]["qg"], [0.0451820, 0.0290373, 0.0343748]; atol=1e-3))
        # end

        @testset "5-bus coupled meshed infeasible acp opf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_b.m")
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_INFEASIBLE
        end

        @testset "5-bus coupled radial no shunt acp opf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_c_r_a.m")
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 55436.1; atol=1e-1)

            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c], 0.4; atol=1e-3) for c in 1:mp_data["conductors"])
            @test all(isapprox(result["solution"]["bus"]["2"]["vm"][c], 1.08564; atol=1e-3) for c in 1:mp_data["conductors"])
            @test all(isapprox(result["solution"]["bus"]["2"]["va"][c], PMD._wrap_to_pi(0.04905-2*pi/mp_data["conductors"]*(c-1)); atol=1e-3) for c in 1:mp_data["conductors"])
        end

        @testset "5-bus coupled radial shunt acp opf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_c_r_b.m")
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 56075.1; atol=1e-1)

            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.4; atol=1e-3) for c in 1:mp_data["conductors"])
            @test all(isapprox(result["solution"]["bus"]["2"]["va"][c], PMD._wrap_to_pi(0.055338-2*pi/mp_data["conductors"]*(c-1)); atol=5e-3) for c in 1:mp_data["conductors"])
        end

        # @testset "voltage unbalance constrained acp opf" begin
        #     pmd_data = PMD.parse_file("../test/data/matlab/case_bctr.m")
        #     # We check the equations by comparing against the value calculated by the solution
        #     # builder for the active constraint
        #     constr_keys = ["vm_vuf_max", "vm_seq_neg_max", "vm_seq_zero_max", "vm_seq_pos_max", "vm_ll_max", "vm_ll_min"]
        #     constr_lims = [0.04, 0.04, 0.04, 1.02, PMD.MultiConductorVector(ones(3)*1.07), PMD.MultiConductorVector(ones(3)*1.01)]
        #     sol_keys = ["vuf", "vm_seq_neg", "vm_seq_zero", "vm_seq_pos", "vm_ll", "vm_ll"]
        #     for i in 1:length(constr_keys)
        #         pmd = deepcopy(pmd_data)
        #         pmd["bus"]["3"][constr_keys[i]] = constr_lims[i]
        #         sol = PMD.run_mc_opf_bctr(pmd, PMs.ACPPowerModel, ipopt_solver, multiconductor=true)
        #         # the minimum is needed for the LL constraints; only one out of three will be active
        #         @test minimum(abs.(sol["solution"]["bus"]["3"][sol_keys[i]]-constr_lims[i])) <= 1E-5
        #     end
        # end
    end

    @testset "test dropped phases opf" begin
        @testset "4-bus phase drop acp opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case4_phase_drop.dss")
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0182595; atol=1e-4)

            @test all(isapprox.(result["solution"]["gen"]["1"]["pg"], [5.06513e-5, 6.0865e-5, 7.1119e-5]; atol=1e-7))
            @test all(isapprox.(result["solution"]["bus"]["2"]["vm"], [0.990023, 1.000000, 1.000000]; atol=1e-4))
        end

        @testset "4-bus phase drop acr opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case4_phase_drop.dss")
            result = run_mc_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0182595; atol=1e-4)

            calc_vm(id) = sqrt.(result["solution"]["bus"][id]["vr"].^2+result["solution"]["bus"][id]["vi"].^2)
            @test all(isapprox.(result["solution"]["gen"]["1"]["pg"], [5.06513e-5, 6.0865e-5, 7.1119e-5]; atol=1e-7))
            @test isapprox(calc_vm("2")[1], 0.98995; atol=1.5e-4)
        end

        @testset "5-bus phase drop acp opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss")
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0599389; atol=1e-4)

            @test all(isapprox.(result["solution"]["gen"]["1"]["pg"], [0.00015236280779412599, 0.00019836795302238667, 0.0002486642034746673]; atol=1e-7))
            @test all(isapprox.(result["solution"]["bus"]["2"]["vm"], [0.97351, 0.96490, 0.95646]; atol=1e-4))
        end

        @testset "5-bus phase drop acr opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss")
            result = run_mc_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0599400; atol = 1e-4)

            calc_vm(id) = sqrt.(result["solution"]["bus"][id]["vr"].^2+result["solution"]["bus"][id]["vi"].^2)
            @test all(isapprox.(result["solution"]["gen"]["1"]["pg"], [0.00015236280779412599, 0.00019836795302238667, 0.0002486688793741932]; atol=1e-7))
            @test all(isapprox.(calc_vm("2"), [0.9735188343958152, 0.9649003198689144, 0.9564593296045091]; atol=1e-4))
        end

        @testset "5-bus phase drop dcp opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss")
            result = run_mc_opf(mp_data, PMs.DCPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0540021; atol=1e-4)
        end

        @testset "5-bus phase drop nfa opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss")
            result = run_mc_opf(mp_data, PMs.NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.054; atol=1e-4)
        end
    end

    @testset "test opendss opf" begin
        @testset "2-bus diagonal acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case2_diag.dss")
            sol = PMD.run_mc_opf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            @test all(isapprox.(sol["solution"]["bus"]["2"]["vm"], 0.984377; atol=1e-4))
            @test all(isapprox.(sol["solution"]["bus"]["2"]["va"], PMD._wrap_to_pi.([2 * pi / pmd["conductors"] * (1 - c) - deg2rad(0.79) for c in 1:pmd["conductors"]]); atol=deg2rad(0.2)))

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
        end

        @testset "3-bus balanced acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
            sol = PMD.run_mc_opf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            for (bus, va, vm) in zip(["1", "2", "3"], [0.0, deg2rad(-0.03), deg2rad(-0.07)], [0.9959, 0.986973, 0.976605])
                @test all(isapprox.(sol["solution"]["bus"][bus]["va"], PMD._wrap_to_pi.([2 * pi / pmd["conductors"] * (1 - c) + va for c in 1:pmd["conductors"]]); atol=deg2rad(0.01)))
                @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=1e-4))
            end

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018345; atol=1e-6)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00919404; atol=1.2e-5)
        end

        # @testset "3-bus balanced soc opf" begin
        #     pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
        #     sol = PMD.run_mc_opf(pmd, PMs.SOCWRPowerModel, ipopt_solver)
        #
        #     @test sol["termination_status"] == PMs.LOCALLY_SOLVED
        # end

        @testset "3-bus unbalanced acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced.dss")
            sol = PMD.run_mc_opf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            for (bus, va, vm) in zip(["1", "2", "3"],
                                    [0.0, deg2rad.([-0.22, -0.11, 0.12]), deg2rad.([-0.48, -0.24, 0.27])],
                                    [0.9959, [0.980937, 0.98936, 0.987039], [0.963546, 0.981757, 0.976779]])
                @test all(isapprox.(sol["solution"]["bus"][bus]["va"], PMD._wrap_to_pi.([2 * pi / pmd["conductors"] * (1 - c) for c in 1:pmd["conductors"]]) .+ va; atol=deg2rad(0.01)))
                @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=1e-5))
            end

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214812; atol=1e-6)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00927263; atol=1e-5)
        end

        # @testset "3-bus unbalanced soc opf" begin
        #     pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced.dss")
        #     sol = PMD.run_mc_opf(pmd, PMs.SOCWRPowerModel, ipopt_solver)
        #
        #     @test sol["termination_status"] == PMs.LOCALLY_SOLVED
        # end

        @testset "3-bus unbalanced isc acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced_isc.dss")
            sol = PMD.run_mc_opf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.0183961; atol=1e-4)
        end

        @testset "3-bus balanced pv acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced_pv.dss")

            @test length(pmd["gen"]) == 2
            @test all(pmd["gen"]["2"]["qmin"] .== -pmd["gen"]["2"]["qmax"])
            @test all(pmd["gen"]["2"]["pmax"] .== pmd["gen"]["2"]["qmax"])
            @test all(pmd["gen"]["2"]["pmin"] .== 0.0)

            sol = PMD.run_mc_opf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED
            @test sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]) < 0.0
            @test sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]) < 0.0
            @test isapprox(sum(sol["solution"]["gen"]["2"]["pg"] * sol["solution"]["baseMVA"]), 0.0182074; atol=1e-4)
            @test isapprox(sum(sol["solution"]["gen"]["2"]["qg"] * sol["solution"]["baseMVA"]), 0.00919404; atol=1e-4)
        end

        @testset "3-bus unbalanced single-phase pv acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced_1phase-pv.dss")
            sol = PMD.run_mc_opf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0196116; atol=1e-3)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923107; atol=1e-3)

            @test all(sol["solution"]["gen"]["2"]["pg"][2:3] .== 0.0)
            @test all(sol["solution"]["gen"]["2"]["qg"][2:3] .== 0.0)
        end

        @testset "3-bus balanced capacitor acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced_cap.dss")
            sol = PMD.run_mc_pf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            @test all(abs(sol["solution"]["bus"]["3"]["vm"][c]-0.98588)<=1E-4 for c in 1:3)
            @test all(abs(sol["solution"]["bus"]["2"]["vm"][c]-0.99127)<=1E-4 for c in 1:3)
        end

        @testset "3w transformer nfa opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/ut_trans_3w_dyy_1.dss")
            result = run_mc_opf(mp_data, PMs.NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.616; atol=1e-3)
        end

        # @testset "1-bus soc opf" begin
        #     pmd = PMD.parse_file("../test/data/opendss/test_simple.dss")
        #     sol = PMD.run_mc_opf(pmd, PowerModels.SOCWRPowerModel, ipopt_solver)
        #
        #     @test sol["termination_status"] == PMs.LOCALLY_SOLVED
        # end
    end
end
