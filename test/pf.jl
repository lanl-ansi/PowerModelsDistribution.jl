@info "running power flow (pf) tests"

@testset "test pf" begin
    case2_diag = PMD.parse_file("../test/data/opendss/case2_diag.dss"; data_model="mathematical")
    case3_balanced = PMD.parse_file("../test/data/opendss/case3_balanced.dss"; data_model="mathematical")
    case3_unbalanced = PMD.parse_file("../test/data/opendss/case3_unbalanced.dss"; data_model="mathematical")
    case5_phase_drop = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss"; data_model="mathematical")
    case_mxshunt = PMD.parse_file("../test/data/opendss/case_mxshunt.dss"; data_model="mathematical")

    @testset "2-bus diagonal acp pf" begin
        sol = PMD.run_mc_pf(case2_diag, PMs.ACPPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        @test all(isapprox.(sol["solution"]["bus"]["1"]["vm"], 0.984377; atol=1e-4))
        @test all(isapprox.(sol["solution"]["bus"]["1"]["va"], [PMD._wrap_to_pi(2*pi/case2_diag["conductors"]*(1-c) - deg2rad(0.79)) for c in 1:case2_diag["conductors"]]; atol=deg2rad(0.2)))

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
    end

    @testset "2-bus diagonal acr pf" begin
        sol = PMD.run_mc_pf(case2_diag, PMs.ACRPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        calc_vm(id) = sqrt.(sol["solution"]["bus"][id]["vr"].^2+sol["solution"]["bus"][id]["vi"].^2)
        calc_va(id) = atan.(sol["solution"]["bus"][id]["vi"], sol["solution"]["bus"][id]["vr"])
        @test all(isapprox.(calc_vm("1"), 0.984377; atol=1e-4))
        @test all(isapprox.(calc_va("1"), [PMD._wrap_to_pi(2*pi/case2_diag["conductors"]*(1-c) - deg2rad(0.79)) for c in 1:case2_diag["conductors"]]; atol=deg2rad(0.2)))

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
    end

    @testset "3-bus balanced acp pf" begin
        sol = PMD.run_mc_pf(case3_balanced, PMs.ACPPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        for (bus, va, vm) in zip(["2", "1", "3"], [0.0, deg2rad(-0.08), deg2rad(-0.17)], [0.9959, 0.986976, 0.976611])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"], [PMD._wrap_to_pi(2*pi/case3_balanced["conductors"]*(1-c) .+ va) for c in 1:case3_balanced["conductors"]]; atol=deg2rad(0.2)))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=1e-3))
        end

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183456; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923328; atol=1e-4)
    end

    @testset "3-bus balanced acr pf" begin
        sol = PMD.run_mc_pf(case3_balanced, PMs.ACRPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        for (bus, va, vm) in zip(["2", "1", "3"], [0.0, deg2rad(-0.08), deg2rad(-0.17)], [0.9959, 0.986976, 0.976611])
            calc_vm(id) = sqrt.(sol["solution"]["bus"][id]["vr"].^2+sol["solution"]["bus"][id]["vi"].^2)
            calc_va(id) = atan.(sol["solution"]["bus"][id]["vi"], sol["solution"]["bus"][id]["vr"])
            @test all(isapprox.(calc_va(bus), [PMD._wrap_to_pi(2*pi/case3_balanced["conductors"]*(1-c) .+ va) for c in 1:case3_balanced["conductors"]]; atol=deg2rad(0.2)))
            @test all(isapprox.(calc_vm(bus), vm; atol=1e-3))
        end

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183456; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923328; atol=1e-4)
    end

    @testset "3-bus balanced no linecode basefreq defined acp pf" begin
        sol = PMD.run_mc_pf(case3_balanced, PMs.ACPPowerModel, ipopt_solver)

        pmd2 = PMD.parse_file("../test/data/opendss/case3_balanced_basefreq.dss"; data_model="mathematical")
        sol2 = PMD.run_mc_pf(pmd2, PMs.ACPPowerModel, ipopt_solver)

        @test all(all(isapprox.(bus["vm"], sol2["solution"]["bus"][i]["vm"]; atol=1e-8)) for (i, bus) in sol["solution"]["bus"])
        @test all(all(isapprox.(bus["va"], sol2["solution"]["bus"][i]["va"]; atol=1e-8)) for (i, bus) in sol["solution"]["bus"])
        @test all(isapprox(sum(sol["solution"]["gen"]["1"][field] * sol["solution"]["baseMVA"]), sum(sol2["solution"]["gen"]["1"][field] * sol2["solution"]["baseMVA"]); atol=1e-8) for field in ["pg", "qg"])
    end

    @testset "3-bus balanced w/ switch acp pf" begin
        pmd = PMD.parse_file("../test/data/opendss/case3_balanced_switch.dss"; data_model="mathematical")
        sol = PMD.run_mc_pf(pmd, PMs.ACPPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        for (bus, va, vm) in zip(["2", "1", "3"], [0.0, 0.0, deg2rad(-0.04)], [0.9959, 0.995729, 0.985454])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"], [PMD._wrap_to_pi(2*pi/pmd["conductors"]*(1-c) .+ va) for c in 1:pmd["conductors"]]; atol=deg2rad(0.2)))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=1e-3))
        end
    end

    @testset "3-bus unbalanced acp pf" begin
        sol = PMD.run_mc_pf(case3_unbalanced, PMs.ACPPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        for (bus, va, vm) in zip(["2", "1", "3"],
                                [0.0, deg2rad.([-0.22, -0.11, 0.12]), deg2rad.([-0.48, -0.24, 0.27])],
                                [0.9959, [0.98094, 0.989365, 0.987043], [0.96355, 0.981767, 0.976786]])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"], PMD._wrap_to_pi([2*pi/case3_unbalanced["conductors"]*(1-c) for c in 1:case3_unbalanced["conductors"]] .+ va); atol=deg2rad(0.2)))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=2e-3))
        end

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214835; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00932693; atol=1e-4)
    end

    @testset "3-bus unbalanced acr pf" begin
        sol = PMD.run_mc_pf(case3_unbalanced, PMs.ACRPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        for (bus, va, vm) in zip(["2", "1", "3"],
                                [0.0, deg2rad.([-0.22, -0.11, 0.12]), deg2rad.([-0.48, -0.24, 0.27])],
                                [0.9959, [0.98094, 0.989365, 0.987043], [0.96355, 0.981767, 0.976786]])
            calc_vm(id) = sqrt.(sol["solution"]["bus"][id]["vr"].^2+sol["solution"]["bus"][id]["vi"].^2)
            calc_va(id) = atan.(sol["solution"]["bus"][id]["vi"], sol["solution"]["bus"][id]["vr"])
            @test all(isapprox.(calc_va(bus), PMD._wrap_to_pi([2*pi/case3_unbalanced["conductors"]*(1-c) for c in 1:case3_unbalanced["conductors"]] .+ va); atol=deg2rad(0.2)))
            @test all(isapprox.(calc_vm(bus), vm; atol=2e-3))
        end

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214835; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00932693; atol=1e-4)
    end

    @testset "3-bus unbalanced w/ assymetric linecode & phase order swap acp pf" begin
        pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced_assym_swap.dss"; data_model="mathematical")
        sol = PMD.run_ac_mc_pf(pmd, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        @test all(isapprox.(sol["solution"]["bus"]["1"]["vm"], [0.983453, 0.98718, 0.981602]; atol=1e-5))
        @test all(isapprox.(sol["solution"]["bus"]["1"]["va"], deg2rad.([-0.07, -120.19, 120.29]); atol=1e-2))
    end

    @testset "5-bus phase drop acp pf" begin
        result = run_mc_pf(case5_phase_drop, PMs.ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol = 1e-4)

        @test all(isapprox.(result["solution"]["bus"]["3"]["vm"], [0.973519, 0.964902, 0.956465]; atol = 1e-3))
    end

    @testset "5-bus phase drop acr pf" begin
        result = run_mc_pf(case5_phase_drop, PMs.ACRPowerModel, ipopt_solver)

        @test result["termination_status"] == PMs.LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0; atol = 1e-4)

        calc_vm(id) = sqrt.(result["solution"]["bus"][id]["vr"].^2+result["solution"]["bus"][id]["vi"].^2)
        @test isapprox(calc_vm("3")[1], 0.973519; atol = 1e-4)
        @test isapprox(calc_vm("3")[2], 0.964902; atol = 1e-4)
        @test isapprox(calc_vm("3")[3], 0.956465; atol = 1e-4)
    end

    @testset "matrix branch shunts acp pf" begin
        sol = PMD.run_ac_mc_pf(case_mxshunt, ipopt_solver)

        @test all(isapprox.(sol["solution"]["bus"]["2"]["vm"], [0.987399, 0.981300, 1.003536]; atol=1E-6))
    end

    @testset "matrix branch shunts acr pf" begin
        pm = PMs.instantiate_model(case_mxshunt, PMs.ACRPowerModel, PMD.build_mc_pf, ref_extensions=[PMD.ref_add_arcs_trans!], multiconductor=true)
        sol = PMs.optimize_model!(pm, optimizer=ipopt_solver)

        calc_vm(id) = sqrt.(sol["solution"]["bus"][id]["vr"].^2+sol["solution"]["bus"][id]["vi"].^2)
        @test all(isapprox.(calc_vm("2"), [0.987399, 0.981299, 1.003537]; atol=1E-6))
    end

    @testset "virtual sourcebus creation acp pf" begin
        pmd = parse_file("../test/data/opendss/virtual_sourcebus.dss"; data_model="mathematical")
        result = run_ac_mc_pf(pmd, ipopt_solver)
        @test result["termination_status"] == PMs.LOCALLY_SOLVED

        @test all(all(isapprox.(result["solution"]["bus"]["$n"]["vm"], [0.961352, 0.999418, 1.00113]; atol=1e-6)) for n in [1, 2])
        @test all(all(isapprox.(rad2deg.(result["solution"]["bus"]["$n"]["va"]), [-1.25, -120.06, 120.0]; atol=1e-1)) for n in [1, 2])
    end

    @testset "2-bus diagonal ivr pf" begin
        sol = PMD.run_mc_pf_iv(case2_diag, PMs.IVRPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.018209; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.000208979; atol=1e-5)
    end

    @testset "3-bus balanced ivr pf" begin
        sol = PMD.run_mc_pf_iv(case3_balanced, PMs.IVRPowerModel, ipopt_solver)

        @test sol["termination_status"] == PMs.LOCALLY_SOLVED

        @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183456; atol=1e-5)
        @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923328; atol=1e-4)
    end

end
