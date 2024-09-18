@info "running power flow (pf) tests"

@testset "test pf" begin
    @testset "2-bus diagonal acp pf" begin
        result = solve_mc_pf(case2_diag, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"], 0.227339; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["primary"]["va"], [-0.657496, -120.657, 119.343]; atol=1e-2))

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 18.20887; atol=1e-4)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]),  0.20887; atol=1e-4)
    end

    @testset "2-bus diagonal acr pf" begin
        result = solve_mc_pf(case2_diag, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"], 0.227339; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["primary"]["va"], [-0.657496, -120.657, 119.343]; atol=1e-2))

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 18.20887; atol=1e-4)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]),  0.20887; atol=1e-4)
    end

    @testset "2-bus diagonal ivr pf" begin
        result = solve_mc_pf(case2_diag, IVRUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 18.20896; atol=1e-4)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]),  0.20896; atol=1e-4)
    end

    @testset "3-bus balanced acp pf" begin
        result = solve_mc_pf(case3_balanced, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        @test all(isapprox.(result["solution"]["bus"]["sourcebus"]["vm"], 0.229993; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["sourcebus"]["va"], [0.0, -120.0, 120.0]; atol=1e-2))

        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"], 0.227932; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["primary"]["va"], [-0.08, -120.08, 119.92]; atol=0.2))

        @test all(isapprox.(result["solution"]["bus"]["loadbus"]["vm"], 0.225537; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["loadbus"]["va"], [-0.17, -120.17, 119.83]; atol=0.2))

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 18.34478; atol=1e-2)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]),  9.19392; atol=1e-2)
    end

    @testset "3-bus balanced acr pf" begin
        result = solve_mc_pf(case3_balanced, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        @test all(isapprox.(result["solution"]["bus"]["sourcebus"]["vm"], 0.229993; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["sourcebus"]["va"], [0.0, -120.0, 120.0]; atol=1e-2))

        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"], 0.227932; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["primary"]["va"], [-0.08, -120.08, 119.92]; atol=0.2))

        @test all(isapprox.(result["solution"]["bus"]["loadbus"]["vm"], 0.225537; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["loadbus"]["va"], [-0.17, -120.17, 119.83]; atol=0.2))

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 18.34478; atol=1e-2)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]),  9.19392; atol=1e-2)
    end

    @testset "3-bus balanced ivr pf" begin
        result = solve_mc_pf(case3_balanced, IVRUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 18.34498; atol=1e-5)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]),  9.19404; atol=1e-4)
    end

    @testset "3-bus balanced no linecode basefreq defined acp pf" begin
        result1 = solve_mc_pf(case3_balanced, ACPUPowerModel, ipopt_solver)
        result2 = solve_mc_pf(case3_balanced_basefreq, ACPUPowerModel, ipopt_solver)

        @test all(all(isapprox.(result1["solution"]["bus"]["sourcebus"]["vm"], result2["solution"]["bus"]["sourcebus"]["vm"]; atol=1e-8)) for (i, bus) in result1["solution"]["bus"])
        @test all(all(isapprox.(result1["solution"]["bus"]["sourcebus"]["va"], result2["solution"]["bus"]["sourcebus"]["va"]; atol=1e-8)) for (i, bus) in result1["solution"]["bus"])

        @test all(all(isapprox.(result1["solution"]["bus"]["primary"]["vm"], result2["solution"]["bus"]["primary"]["vm"]; atol=1e-8)) for (i, bus) in result1["solution"]["bus"])
        @test all(all(isapprox.(result1["solution"]["bus"]["primary"]["va"], result2["solution"]["bus"]["primary"]["va"]; atol=1e-8)) for (i, bus) in result1["solution"]["bus"])

        @test all(all(isapprox.(result1["solution"]["bus"]["loadbus"]["vm"], result2["solution"]["bus"]["loadbus"]["vm"]; atol=1e-8)) for (i, bus) in result1["solution"]["bus"])
        @test all(all(isapprox.(result1["solution"]["bus"]["loadbus"]["va"], result2["solution"]["bus"]["loadbus"]["va"]; atol=1e-8)) for (i, bus) in result1["solution"]["bus"])

        @test isapprox(sum(result1["solution"]["voltage_source"]["source"]["pg"]), sum(result2["solution"]["voltage_source"]["source"]["pg"]); atol=1e-8)
        @test isapprox(sum(result1["solution"]["voltage_source"]["source"]["qg"]), sum(result2["solution"]["voltage_source"]["source"]["qg"]); atol=1e-8)
    end

    @testset "3-bus unbalanced acp pf" begin
        result = solve_mc_pf(case3_unbalanced, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case3_unbalanced["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["sourcebus"]["vm"] ./ vbase, [0.9959, 0.9959, 0.9959]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["sourcebus"]["va"], [0.0, -120.0, 120.0]; atol=1e-2))

        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, [0.98094, 0.989365, 0.987043]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["primary"]["va"], [-0.22, -120.11, 120.12]; atol=1e-2))

        @test all(isapprox.(result["solution"]["bus"]["loadbus"]["vm"] ./ vbase, [0.96355, 0.981767, 0.976786]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["loadbus"]["va"], [-0.48, -120.24, 120.27]; atol=1e-2))

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 21.4812; atol=1e-2)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), 9.27263; atol=1e-2)
    end

    @testset "3-bus unbalanced acr pf" begin
        result = solve_mc_pf(case3_unbalanced, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case3_unbalanced["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["sourcebus"]["vm"] ./ vbase, [0.9959, 0.9959, 0.9959]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["sourcebus"]["va"], [0.0, -120.0, 120.0]; atol=1e-2))

        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, [0.98094, 0.989365, 0.987043]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["primary"]["va"], [-0.22, -120.11, 120.12]; atol=1e-2))

        @test all(isapprox.(result["solution"]["bus"]["loadbus"]["vm"] ./ vbase, [0.96355, 0.981767, 0.976786]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["loadbus"]["va"], [-0.48, -120.24, 120.27]; atol=1e-2))

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 21.4812; atol=1e-2)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), 9.27263; atol=1e-2)
    end

    @testset "3-bus unbalanced w/ asymmetric linecode & phase order swap acp pf" begin
        result = solve_mc_pf(case3_unbalanced_assym_swap, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case3_unbalanced_assym_swap["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["primary"]["vm"] ./ vbase, [0.983453, 0.98718, 0.981602]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"]["primary"]["va"], [-0.07, -120.19, 120.29]; atol=1e-2))
    end

    @testset "5-bus phase drop acp pf" begin
        result = solve_mc_pf(case5_phase_drop, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case5_phase_drop["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["midbus"]["vm"] ./ vbase, [0.973519, 0.964902, 0.956465]; atol = 1e-4))
    end

    @testset "5-bus phase drop acr pf" begin
        result = solve_mc_pf(case5_phase_drop, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        vbase = case5_phase_drop["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["midbus"]["vm"] ./ vbase, [0.973519, 0.964902, 0.956465]; atol=1e-4))
    end

    @testset "matrix branch shunts acp pf" begin
        result = solve_mc_pf(case2_mxshunt, ACPUPowerModel, ipopt_solver)

        vbase = case2_mxshunt["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["loadbus"]["vm"] ./ vbase, [0.987399, 0.981300, 1.003536]; atol=1e-4))
    end

    @testset "matrix branch shunts acr pf" begin
        result = solve_mc_pf(case2_mxshunt, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        vbase = case2_mxshunt["settings"]["vbases_default"]["sourcebus"]
        @test all(isapprox.(result["solution"]["bus"]["loadbus"]["vm"] ./ vbase, [0.987399, 0.981300, 1.003536]; atol=1e-4))
    end

    @testset "virtual sourcebus creation acp pf" begin
        result = solve_mc_pf(transform_data_model(case2_virtual_sourcebus), ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        @test all(all(isapprox.(result["solution"]["bus"]["$n"]["vm"], [0.961352, 0.999418, 1.00113]; atol=1e-4)) for n in [1, 2])
        @test all(all(isapprox.(result["solution"]["bus"]["$n"]["va"], deg2rad.([-1.25, -120.06, 120.0]); atol=1e-2)) for n in [1, 2])
    end
end
