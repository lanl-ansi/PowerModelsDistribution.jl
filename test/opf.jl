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
            @test all(isapprox(result["solution"]["bus"]["2"]["va"][c], -0.0538204; atol=1e-5) for c in 1:mp_data["conductors"])
        end

        @testset "5-bus matpower acr opf" begin
            mp_data = PMs.parse_file("../test/data/matpower/case5.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 45522.096; atol=1e-1)

            calc_va(id) = atan.(result["solution"]["bus"][id]["vi"], result["solution"]["bus"][id]["vr"])
            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.3999999; atol=1e-3) for c in 1:mp_data["conductors"])
            @test all(isapprox(calc_va("2")[c], -0.0538204; atol=1e-5) for c in 1:mp_data["conductors"])
        end

        @testset "30-bus matpower acp opf" begin
            mp_data = PMs.parse_file("../test/data/matpower/case30.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 614.007; atol=1e-1)

            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  2.192189; atol=1e-3) for c in 1:mp_data["conductors"])
            @test all(isapprox(result["solution"]["bus"]["2"]["va"][c], -0.071853; atol=1e-4) for c in 1:mp_data["conductors"])
        end

        @testset "30-bus matpower acr opf" begin
            mp_data = PMs.parse_file("../test/data/matpower/case30.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 614.007; atol=1e-1)

            calc_va(id) = atan.(result["solution"]["bus"][id]["vi"], result["solution"]["bus"][id]["vr"])
            @test all(isapprox(result["solution"]["gen"]["1"]["pg"][c],  2.192189; atol=1e-3) for c in 1:mp_data["conductors"])
            @test all(isapprox(calc_va("2")[c], -0.071853; atol=1e-4) for c in 1:mp_data["conductors"])
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
    end

    @testset "test dropped phases opf" begin
        @testset "4-bus phase drop acp opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case4_phase_drop.dss")
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0182595; atol=1e-4)

            @test all(isapprox.(result["solution"]["gen"]["1"]["pg"], [5.06513e-5, 6.0865e-5, 7.1119e-5]; atol=1e-7))
            @test isapprox(result["solution"]["bus"]["4"]["vm"][1], 0.98995; atol=1.5e-4)
        end

        @testset "4-bus phase drop acr opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case4_phase_drop.dss")
            result = run_mc_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0182595; atol=1e-4)

            calc_vm(id) = sqrt.(result["solution"]["bus"][id]["vr"].^2+result["solution"]["bus"][id]["vi"].^2)
            @test all(isapprox.(result["solution"]["gen"]["1"]["pg"], [5.06513e-5, 6.0865e-5, 7.1119e-5]; atol=1e-7))
            @test isapprox(calc_vm("4")[1], 0.98995; atol=1.5e-4)
        end

        @testset "5-bus phase drop acp opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss")
            result = run_mc_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0599389; atol=1e-4)

            @test all(isapprox.(result["solution"]["gen"]["1"]["pg"], [0.00015236280779412599, 0.00019836795302238667, 0.0002486642034746673]; atol=1e-7))
            @test all(isapprox.(result["solution"]["bus"]["3"]["vm"], [0.97351, 0.96490, 0.95646]; atol=1e-4))
        end

        @testset "5-bus phase drop acr opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss")
            result = run_mc_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0599400; atol = 1e-4)

            calc_vm(id) = sqrt.(result["solution"]["bus"][id]["vr"].^2+result["solution"]["bus"][id]["vi"].^2)
            @test all(isapprox.(result["solution"]["gen"]["1"]["pg"], [0.00015236280779412599, 0.00019836795302238667, 0.0002486688793741932]; atol=1e-7))
            @test all(isapprox.(calc_vm("3"), [0.9735188343958152, 0.9649003198689144, 0.9564593296045091]; atol=1e-4))
        end

        @testset "5-bus phase drop dcp opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss")
            result = run_mc_opf(mp_data, PMs.DCPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0544220; atol=1e-4)
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

            @test all(isapprox.(sol["solution"]["bus"]["1"]["vm"], 0.984377; atol=1e-4))
            @test all(isapprox.(sol["solution"]["bus"]["1"]["va"], PMD._wrap_to_pi.([2 * pi / pmd["conductors"] * (1 - c) - deg2rad(0.79) for c in 1:pmd["conductors"]]); atol=deg2rad(0.2)))

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0181409; atol=1e-5)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.0; atol=1e-4)
        end

        @testset "3-bus balanced acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
            sol = PMD.run_mc_opf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            for (bus, va, vm) in zip(["2", "1", "3"], [0.0, deg2rad(-0.03), deg2rad(-0.07)], [0.9959, 0.986973, 0.976605])
                @test all(isapprox.(sol["solution"]["bus"][bus]["va"], PMD._wrap_to_pi.([2 * pi / pmd["conductors"] * (1 - c) + va for c in 1:pmd["conductors"]]); atol=deg2rad(0.01)))
                @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=1e-4))
            end

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018276; atol=1e-6)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.008922; atol=1.2e-5)
        end

        @testset "3-bus unbalanced acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced.dss")
            sol = PMD.run_mc_opf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            for (bus, va, vm) in zip(["2", "1", "3"],
                                    [0.0, deg2rad.([-0.22, -0.11, 0.12]), deg2rad.([-0.48, -0.24, 0.27])],
                                    [0.9959, [0.980937, 0.98936, 0.987039], [0.963546, 0.981757, 0.976779]])
                @test all(isapprox.(sol["solution"]["bus"][bus]["va"], PMD._wrap_to_pi.([2 * pi / pmd["conductors"] * (1 - c) for c in 1:pmd["conductors"]]) .+ va; atol=deg2rad(0.01)))
                @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=1e-5))
            end

            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214812; atol=1e-6)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00927263; atol=1e-5)
        end

        @testset "3-bus unbalanced isc acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced_isc.dss")
            sol = PMD.run_mc_opf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.0185; atol=1e-4)
        end

        @testset "3-bus balanced pv acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced_pv.dss")

            @test length(pmd["gen"]) == 2
            @test all(pmd["gen"]["1"]["qmin"] .== -pmd["gen"]["1"]["qmax"])
            @test all(pmd["gen"]["1"]["pmax"] .==  pmd["gen"]["1"]["qmax"])
            @test all(pmd["gen"]["1"]["pmin"] .== 0.0)

            sol = PMD.run_mc_opf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED
            @test sum(sol["solution"]["gen"]["2"]["pg"] * sol["solution"]["baseMVA"]) < 0.0
            @test sum(sol["solution"]["gen"]["2"]["qg"] * sol["solution"]["baseMVA"]) < 0.0
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183685; atol=1e-4)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00919404; atol=1e-4)
        end

        @testset "3-bus unbalanced single-phase pv acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced_1phase-pv.dss")
            sol = PMD.run_mc_opf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            @test isapprox(sum(sol["solution"]["gen"]["2"]["pg"] * sol["solution"]["baseMVA"]), 0.01838728; atol=1e-3)
            @test isapprox(sum(sol["solution"]["gen"]["2"]["qg"] * sol["solution"]["baseMVA"]), 0.00756634; atol=1e-3)

            @test all(sol["solution"]["gen"]["1"]["pg"][2:3] .== 0.0)
            @test all(sol["solution"]["gen"]["1"]["qg"][2:3] .== 0.0)
        end

        @testset "3-bus balanced capacitor acp opf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced_cap.dss")
            sol = PMD.run_mc_pf(pmd, PMs.ACPPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED

            @test all(abs(sol["solution"]["bus"]["3"]["vm"][c]-0.98588)<=1E-4 for c in 1:3)
            @test all(abs(sol["solution"]["bus"]["1"]["vm"][c]-0.99127)<=1E-4 for c in 1:3)
        end

        @testset "3w transformer nfa opf" begin
            mp_data = PMD.parse_file("../test/data/opendss/ut_trans_3w_dyy_1.dss")
            result = run_mc_opf(mp_data, PMs.NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.616; atol=1e-3)
        end
    end
end
