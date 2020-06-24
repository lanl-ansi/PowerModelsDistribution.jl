@info "running optimal power flow (opf) tests"

@testset "test opf" begin
    @testset "test matpower opf" begin
        case5 = PM.parse_file("../test/data/matpower/case5.m")
        case30 = PM.parse_file("../test/data/matpower/case30.m")

        make_multiconductor!(case5, 3)
        make_multiconductor!(case30, 3)

        @testset "5-bus matpower acp opf" begin
            result = run_mc_opf(case5, ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 45522.096; atol=1e-1)

            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.3999999; atol=1e-3) for c in 1:case5["conductors"])
            @test all(isapprox(result["solution"]["bus"]["2"]["va"][c], -0.0538204; atol=1e-5) for c in 1:case5["conductors"])
        end

        @testset "5-bus matpower acr opf" begin
            result = run_mc_opf(case5, ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 45522.096; atol=1e-1)

            calc_va(id) = atan.(result["solution"]["bus"][id]["vi"], result["solution"]["bus"][id]["vr"])
            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.3999999; atol=1e-3) for c in 1:case5["conductors"])
            @test all(isapprox(calc_va("2")[c], -0.0538204; atol=1e-5) for c in 1:case5["conductors"])
        end

        @testset "30-bus matpower acp opf" begin
            result = run_mc_opf(case30, ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 614.007; atol=1e-1)

            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  2.192189; atol=1e-3) for c in 1:case30["conductors"])
            @test all(isapprox(result["solution"]["bus"]["2"]["va"][c], -0.071853; atol=1e-4) for c in 1:case30["conductors"])
        end

        @testset "30-bus matpower acr opf" begin
            result = run_mc_opf(case30, ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 614.007; atol=1e-1)

            calc_va(id) = atan.(result["solution"]["bus"][id]["vi"], result["solution"]["bus"][id]["vr"])
            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  2.192189; atol=1e-3) for c in 1:case30["conductors"])
            @test all(isapprox(calc_va("2")[c], -0.071853; atol=1e-4) for c in 1:case30["conductors"])
        end

        @testset "30-bus matpower dcp opf" begin
            result = run_mc_opf(case30, DCPPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 566.112; atol=1e-1)
        end

        @testset "30-bus matpower nfa opf" begin
            result = run_mc_opf(case30, NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 458.006; atol=1e-1)
        end
    end

    @testset "test dropped phases opf" begin
        case4_phase_drop = parse_file("../test/data/opendss/case4_phase_drop.dss")
        case5_phase_drop = parse_file("../test/data/opendss/case5_phase_drop.dss")

        @testset "4-bus phase drop acp opf" begin
            result = run_mc_opf(case4_phase_drop, ACPPowerModel, ipopt_solver, make_si=false)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0182595; atol=1e-4)

            @test all(isapprox.(result["solution"]["voltage_source"]["source"]["pg"], [5.06513e-5, 6.0865e-5, 7.1119e-5]; atol=1e-7))
            @test isapprox(result["solution"]["bus"]["loadbus1"]["vm"][1], 0.98995; atol=1.5e-4)
        end

        @testset "4-bus phase drop acr opf" begin
            result = run_mc_opf(case4_phase_drop, ACRPowerModel, ipopt_solver; make_si=false)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0182595; atol=1e-4)

            @test all(isapprox.(result["solution"]["voltage_source"]["source"]["pg"], [5.06513e-5, 6.0865e-5, 7.1119e-5]; atol=1e-7))
            @test isapprox(calc_vm_acr(result, "loadbus1")[1], 0.98995; atol=1.5e-4)
        end

        @testset "5-bus phase drop acp opf" begin
            result = run_mc_opf(case5_phase_drop, ACPPowerModel, ipopt_solver; make_si=false)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0599389; atol=1e-4)

            @test all(isapprox.(result["solution"]["voltage_source"]["source"]["pg"], [0.00015236280779412599, 0.00019836795302238667, 0.0002486642034746673]; atol=1e-7))
            @test all(isapprox.(result["solution"]["bus"]["midbus"]["vm"], [0.97351, 0.96490, 0.95646]; atol=1e-4))
        end

        @testset "5-bus phase drop acr opf" begin
            result = run_mc_opf(case5_phase_drop, ACRPowerModel, ipopt_solver; make_si=false)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0599400; atol = 1e-4)

            @test all(isapprox.(result["solution"]["voltage_source"]["source"]["pg"], [0.00015236280779412599, 0.00019836795302238667, 0.0002486688793741932]; atol=1e-7))
            @test all(isapprox.(calc_vm_acr(result, "midbus"), [0.9735188343958152, 0.9649003198689144, 0.9564593296045091]; atol=1e-4))
        end

        @testset "5-bus phase drop dcp opf" begin
            result = run_mc_opf(case5_phase_drop, DCPPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0544220; atol=1e-4)
        end

        @testset "5-bus phase drop nfa opf" begin
            result = run_mc_opf(case5_phase_drop, NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.054; atol=1e-4)
        end
    end

    @testset "test opendss opf" begin
        @testset "2-bus diagonal acp opf" begin
            pmd = parse_file("../test/data/opendss/case2_diag.dss")
            sol = run_mc_opf(pmd, ACPPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED

            @test all(isapprox.(sol["solution"]["bus"]["primary"]["vm"], 0.984377; atol=1e-4))
            @test all(isapprox.(sol["solution"]["bus"]["primary"]["va"], [0, -120, 120] .- 0.79; atol=0.2))

            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * sol["solution"]["settings"]["sbase"]), 0.0181409; atol=1e-5)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * sol["solution"]["settings"]["sbase"]), 0.0; atol=1e-4)
        end

        @testset "3-bus balanced acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_balanced.dss")
            sol = run_mc_opf(pmd, ACPPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED

            for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"], [0.0, -0.03, -0.07], [0.9959, 0.986973, 0.976605])
                @test all(isapprox.(sol["solution"]["bus"][bus]["va"], [0, -120, 120] .+ va; atol=0.01))
                @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=1e-4))
            end

            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * sol["solution"]["settings"]["sbase"]), 0.018276; atol=1e-6)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * sol["solution"]["settings"]["sbase"]), 0.008922; atol=1.2e-5)
        end

        @testset "3-bus unbalanced acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_unbalanced.dss")
            sol = run_mc_opf(pmd, ACPPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED

            for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"],
                                    [0.0, [-0.22, -0.11, 0.12], [-0.48, -0.24, 0.27]],
                                    [0.9959, [0.980937, 0.98936, 0.987039], [0.963546, 0.981757, 0.976779]])
                @test all(isapprox.(sol["solution"]["bus"][bus]["va"], [0, -120, 120] .+ va; atol=0.01))
                @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=1e-5))
            end

            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * sol["solution"]["settings"]["sbase"]), 0.0214812; atol=1e-6)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * sol["solution"]["settings"]["sbase"]), 0.00927263; atol=1e-5)
        end

        @testset "3-bus unbalanced isc acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_balanced_isc.dss")
            sol = run_mc_opf(pmd, ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.0185; atol=1e-4)
        end

        @testset "3-bus balanced pv acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_balanced_pv.dss")

            sol = run_mc_opf(pmd, ACPPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED
            @test sum(sol["solution"]["voltage_source"]["source"]["pg"] * sol["solution"]["settings"]["sbase"]) < 0.0
            @test sum(sol["solution"]["voltage_source"]["source"]["qg"] * sol["solution"]["settings"]["sbase"]) < 0.005
            @test isapprox(sum(sol["solution"]["solar"]["pv1"]["pg"] * sol["solution"]["settings"]["sbase"]), 0.0183685; atol=1e-4)
            @test isapprox(sum(sol["solution"]["solar"]["pv1"]["qg"] * sol["solution"]["settings"]["sbase"]), 0.0091449; atol=1e-4)
        end

        @testset "3-bus unbalanced single-phase pv acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_unbalanced_1phase-pv.dss")
            sol = run_mc_opf(pmd, ACPPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED

            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * sol["solution"]["settings"]["sbase"]), 0.01838728; atol=1e-3)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * sol["solution"]["settings"]["sbase"]), 0.00756634; atol=1e-3)

            @test all(sol["solution"]["solar"]["pv1"]["pg"][2:3] .== 0.0)
            @test all(sol["solution"]["solar"]["pv1"]["qg"][2:3] .== 0.0)
        end

        @testset "3-bus balanced capacitor acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_balanced_cap.dss")
            sol = run_mc_pf(pmd, ACPPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED

            @test all(abs(sol["solution"]["bus"]["loadbus"]["vm"][c]-0.98588)<=1E-4 for c in 1:3)
            @test all(abs(sol["solution"]["bus"]["primary"]["vm"][c]-0.99127)<=1E-4 for c in 1:3)
        end

        @testset "3w transformer nfa opf" begin
            mp_data = parse_file("../test/data/opendss/ut_trans_3w_dyy_1.dss")
            result = run_mc_opf(mp_data, NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.616; atol=1e-3)

        end

        @testset "vm and va start with ivr and acr opf" begin

            pmd_eng = parse_file("../test/data/opendss/case3_unbalanced.dss")
            pmd_math = transform_data_model(pmd_eng)

            sol_ivr = run_mc_opf(pmd_math, IVRPowerModel, ipopt_solver)
            sol_acr = run_mc_opf(pmd_math, ACRPowerModel, ipopt_solver)

            pmd_math["bus"]["4"]["vm_start"] = pmd_math["bus"]["4"]["vm"]
            pmd_math["bus"]["4"]["va_start"] = pmd_math["bus"]["4"]["va"]
            pmd_math["bus"]["2"]["vm_start"] = [0.9959, 0.9959, 0.9959]
            pmd_math["bus"]["2"]["va_start"] = [0.00, -2.0944, 2.0944]

            sol_ivr_with_start = run_mc_opf(pmd_math, IVRPowerModel, ipopt_solver)
            sol_acr_with_start = run_mc_opf(pmd_math, ACRPowerModel, ipopt_solver)

            @test isapprox(sol_ivr["objective"], sol_ivr_with_start["objective"]; atol=1e-5)
            @test isapprox(sol_acr["objective"], sol_acr_with_start["objective"]; atol=1e-5)
        end
    end
end
