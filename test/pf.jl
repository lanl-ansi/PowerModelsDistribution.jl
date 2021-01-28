@info "running power flow (pf) tests"

@testset "test pf" begin
    case2_diag = parse_file("../test/data/opendss/case2_diag.dss")
    case3_balanced = parse_file("../test/data/opendss/case3_balanced.dss")
    case3_unbalanced = parse_file("../test/data/opendss/case3_unbalanced.dss")
    case5_phase_drop = parse_file("../test/data/opendss/case5_phase_drop.dss")
    case_mxshunt = parse_file("../test/data/opendss/case_mxshunt.dss")

    @testset "2-bus diagonal acp pf" begin
        sol = run_mc_pf(case2_diag, ACPPowerModel, ipopt_solver)

        @test sol["termination_status"] == LOCALLY_SOLVED

        @test all(isapprox.(sol["solution"]["bus"]["primary"]["vm"], 0.227339; atol=1e-4))
        @test all(isapprox.(sol["solution"]["bus"]["primary"]["va"], [-0.657496, -120.657, 119.343]; atol=0.2))

        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 18.20887; atol=1e-5)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]),  0.20887; atol=1e-5)
    end

    @testset "2-bus diagonal acr pf" begin
        sol = run_mc_pf(case2_diag, ACRPowerModel, ipopt_solver)

        @test sol["termination_status"] == LOCALLY_SOLVED

        @test all(isapprox.(calc_vm_acr(sol, "primary"), 0.227339; atol=1e-4))
        @test all(isapprox.(calc_va_acr(sol, "primary"), [-0.657496, -120.657, 119.343]; atol=0.2))

        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 18.20888; atol=1e-5)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]),  0.20888; atol=1e-5)
    end

    @testset "3-bus balanced acp pf" begin
        sol = run_mc_pf(case3_balanced, ACPPowerModel, ipopt_solver)

        @test sol["termination_status"] == LOCALLY_SOLVED

        for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"], [0.0, -0.08, -0.17], [0.229993, 0.227932, 0.225537])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"], [0, -120, 120] .+ va; atol=0.2))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=1e-3))
        end

        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 18.34478; atol=1e-5)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]),  9.19392; atol=1e-4)
    end

    @testset "3-bus balanced acr pf" begin
        sol = run_mc_pf(case3_balanced, ACRPowerModel, ipopt_solver)

        @test sol["termination_status"] == LOCALLY_SOLVED

        for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"], [0.0, -0.08, -0.17], [0.229993, 0.227932, 0.225537])
            @test all(isapprox.(calc_va_acr(sol, bus), [0, -120, 120] .+ va; atol=0.2))
            @test all(isapprox.(calc_vm_acr(sol, bus), vm; atol=1e-3))
        end

        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 18.34478; atol=1e-5)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]),  9.19392; atol=1e-4)
    end

    @testset "3-bus balanced no linecode basefreq defined acp pf" begin
        sol = run_mc_pf(case3_balanced, ACPPowerModel, ipopt_solver)

        pmd2 = parse_file("../test/data/opendss/case3_balanced_basefreq.dss")
        sol2 = run_mc_pf(pmd2, ACPPowerModel, ipopt_solver)

        @test all(all(isapprox.(bus["vm"], sol2["solution"]["bus"][i]["vm"]; atol=1e-8)) for (i, bus) in sol["solution"]["bus"])
        @test all(all(isapprox.(bus["va"], sol2["solution"]["bus"][i]["va"]; atol=1e-8)) for (i, bus) in sol["solution"]["bus"])
        @test all(isapprox(sum(sol["solution"]["voltage_source"]["source"][field]), sum(sol2["solution"]["voltage_source"]["source"][field]); atol=1e-8) for field in ["pg", "qg"])
    end

    @testset "3-bus unbalanced acp pf" begin
        sol = run_mc_pf(case3_unbalanced, ACPPowerModel, ipopt_solver; make_si=false)

        @test sol["termination_status"] == LOCALLY_SOLVED

        for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"],
                                [0.0, [-0.22, -0.11, 0.12], [-0.48, -0.24, 0.27]],
                                [0.9959, [0.98094, 0.989365, 0.987043], [0.96355, 0.981767, 0.976786]])
            @test all(isapprox.(sol["solution"]["bus"][bus]["va"], [0, -120, 120] .+ va; atol=0.2))
            @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=2e-3))
        end

        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 0.04296; atol=1e-5)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]), 0.01854; atol=1e-4)
    end

    @testset "3-bus unbalanced acr pf" begin
        sol = run_mc_pf(case3_unbalanced, ACRPowerModel, ipopt_solver; make_si=false)

        @test sol["termination_status"] == LOCALLY_SOLVED

        for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"],
                                [0.0, [-0.22, -0.11, 0.12], [-0.48, -0.24, 0.27]],
                                [0.9959, [0.98094, 0.989365, 0.987043], [0.96355, 0.981767, 0.976786]])
            @test all(isapprox.(calc_va_acr(sol, bus), [0, -120, 120] .+ va; atol=0.2))
            @test all(isapprox.(calc_vm_acr(sol, bus), vm; atol=2e-3))
        end

        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 0.04296; atol=1e-5)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]), 0.01854; atol=1e-4)
    end

    @testset "3-bus unbalanced w/ asymmetric linecode & phase order swap acp pf" begin
        pmd = parse_file("../test/data/opendss/case3_unbalanced_assym_swap.dss")
        sol = run_ac_mc_pf(pmd, ipopt_solver; make_si=false)

        @test sol["termination_status"] == LOCALLY_SOLVED

        @test all(isapprox.(sol["solution"]["bus"]["primary"]["vm"], [0.983453, 0.98718, 0.981602]; atol=1e-5))
        @test all(isapprox.(sol["solution"]["bus"]["primary"]["va"], [-0.07, -120.19, 120.29]; atol=1e-2))
    end

    @testset "5-bus phase drop acp pf" begin
        result = run_mc_pf(case5_phase_drop, ACPPowerModel, ipopt_solver; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED

        @test all(isapprox.(result["solution"]["bus"]["midbus"]["vm"], [0.973519, 0.964902, 0.956465]; atol = 1e-3))
    end

    @testset "5-bus phase drop acr pf" begin
        sol = run_mc_pf(case5_phase_drop, ACRPowerModel, ipopt_solver; make_si=false)

        @test sol["termination_status"] == LOCALLY_SOLVED

        @test all(isapprox.(calc_vm_acr(sol, "midbus"), [0.973519, 0.964902, 0.956465]; atol=1e-4))
    end

    @testset "matrix branch shunts acp pf" begin
        sol = run_ac_mc_pf(case_mxshunt, ipopt_solver; make_si=false)

        @test all(isapprox.(sol["solution"]["bus"]["loadbus"]["vm"], [0.987399, 0.981300, 1.003536]; atol=1E-6))
    end

    @testset "matrix branch shunts acr pf" begin
        sol = run_mc_pf(case_mxshunt, ACRPowerModel, ipopt_solver; make_si=false)

        @test all(isapprox.(calc_vm_acr(sol, "loadbus"), [0.987399, 0.981299, 1.003537]; atol=1E-6))
    end

    @testset "virtual sourcebus creation acp pf" begin
        pmd = parse_file("../test/data/opendss/virtual_sourcebus.dss"; data_model=MATHEMATICAL)
        result = run_ac_mc_pf(pmd, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        @test all(all(isapprox.(result["solution"]["bus"]["$n"]["vm"], [0.961352, 0.999418, 1.00113]; atol=1e-6)) for n in [1, 2])
        @test all(all(isapprox.(result["solution"]["bus"]["$n"]["va"], deg2rad.([-1.25, -120.06, 120.0]); atol=1e-1)) for n in [1, 2])
    end

    @testset "2-bus diagonal ivr pf" begin
        sol = run_mc_pf(case2_diag, IVRPowerModel, ipopt_solver)

        @test sol["termination_status"] == LOCALLY_SOLVED

        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 18.20896; atol=1e-5)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]),  0.20896; atol=1e-5)
    end

    @testset "3-bus balanced ivr pf" begin
        sol = run_mc_pf(case3_balanced, IVRPowerModel, ipopt_solver)

        @test sol["termination_status"] == LOCALLY_SOLVED

        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 18.34498; atol=1e-5)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]),  9.19404; atol=1e-4)
    end

end
