@info "running capacitor control tests"

@testset "capacitor control" begin
    @testset "capcontrol_acp" begin
        result = solve_mc_opf_capc(IEEE13_CapControl, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 405.556; atol=5)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), -527.15; atol=200)

        vbase,_ = calc_voltage_bases(IEEE13_CapControl, IEEE13_CapControl["settings"]["vbases_default"])
        @test all(isapprox.(result["solution"]["bus"]["646"]["vm"] ./ vbase["646"], [1.03533, 1.06987]; atol=2e-2))
        @test all(isapprox.(result["solution"]["bus"]["646"]["va"], [-90.1768, 148.498]; atol=3e-1))

        @test all(isapprox.(result["solution"]["shunt"]["c1"]["cap_state"], [1.0]; atol=2e-1))
    end

    @testset "capcontrol_acr" begin
        result = solve_mc_opf_capc(IEEE13_CapControl, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == ALMOST_LOCALLY_SOLVED || result["termination_status"] == LOCALLY_SOLVED

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 405.556; atol=5)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), -527.15; atol=300)

        vbase,_ = calc_voltage_bases(IEEE13_CapControl, IEEE13_CapControl["settings"]["vbases_default"])
        @test all(isapprox.(result["solution"]["bus"]["646"]["vm"] ./ vbase["646"], [1.03533, 1.06987]; atol=2e-2))
        @test all(isapprox.(result["solution"]["bus"]["646"]["va"], [-90.1768, 148.498]; atol=3e-1))

        @test all(isapprox.(result["solution"]["shunt"]["c1"]["cap_state"], [1.0]; atol=2e-1))
    end

    @testset "capcontrol_fbs" begin
        result = solve_mc_opf_capc(IEEE13_CapControl, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 404.784; atol=5)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), -328.146; atol=300)

        vbase,_ = calc_voltage_bases(IEEE13_CapControl, IEEE13_CapControl["settings"]["vbases_default"])
        @test all(isapprox.(result["solution"]["bus"]["646"]["vm"] ./ vbase["646"], [1.03928, 1.05688]; atol=2e-1))
        @test all(isapprox.(result["solution"]["bus"]["646"]["va"], [-89.997, 148.711]; atol=8e-1))

        @test all(isapprox.(result["solution"]["shunt"]["c1"]["cap_state"], [1.0]; atol=6e-1))
    end

    @testset "capcontrol_lpubfdiag" begin
        result = solve_mc_opf_capc(IEEE13_CapControl, LPUBFDiagPowerModel, ipopt_solver_adaptive; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 404.784; atol=5)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), -328.146; atol=900)

        vbase,_ = calc_voltage_bases(IEEE13_CapControl, IEEE13_CapControl["settings"]["vbases_default"])
        @test all(isapprox.(result["solution"]["bus"]["646"]["vm"] ./ vbase["646"], [1.03928, 1.05688]; atol=2e-1))

        @test all(isapprox.(result["solution"]["shunt"]["c1"]["cap_state"], [1.0]; atol=6e-1))
    end

    @testset "capcontrol_fotr" begin
        result = solve_mc_opf_capc(IEEE13_CapControl, FOTRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 404.784; atol=5)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), -328.146; atol=400)

        vbase,_ = calc_voltage_bases(IEEE13_CapControl, IEEE13_CapControl["settings"]["vbases_default"])
        @test all(isapprox.(result["solution"]["bus"]["646"]["vm"] ./ vbase["646"], [1.03928, 1.05688]; atol=2e-1))
        @test all(isapprox.(result["solution"]["bus"]["646"]["va"], [-89.997, 148.711]; atol=1e0))

        @test all(isapprox.(result["solution"]["shunt"]["c1"]["cap_state"], [1.0]; atol=6e-1))
    end
end

