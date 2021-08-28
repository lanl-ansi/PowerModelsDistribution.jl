@info "running capacitor control tests"

@testset "capacitor control" begin
    @testset "capcontrol_acp" begin
        eng = parse_file("../test/data/opendss/IEEE13_CapControl.dss")
        sol = solve_mc_opf_capc(eng, ACPUPowerModel, ipopt_solver; make_si=false)
        @test sol["termination_status"] == LOCALLY_SOLVED
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 0.405556; atol=5e-3)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]), -0.52715; atol=2e-1)
        @test all(isapprox.(sol["solution"]["bus"]["646"]["vm"], [1.03533, 1.06987]; atol=2e-2))
        @test all(isapprox.(sol["solution"]["bus"]["646"]["va"], [-90.1768, 148.498]; atol=3e-1))
        @test all(isapprox.(sol["solution"]["shunt"]["c1"]["cap_state"], [1.0]; atol=2e-1))
    end
    @testset "capcontrol_acr" begin
        eng = parse_file("../test/data/opendss/IEEE13_CapControl.dss")
        sol = solve_mc_opf_capc(eng, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
        @test sol["termination_status"] == LOCALLY_SOLVED
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 0.405556; atol=5e-3)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]), -0.52715; atol=3e-1)
        @test all(isapprox.(sol["solution"]["bus"]["646"]["vm"], [1.03533, 1.06987]; atol=2e-2))
        @test all(isapprox.(sol["solution"]["bus"]["646"]["va"], [-90.1768, 148.498]; atol=3e-1))
        @test all(isapprox.(sol["solution"]["shunt"]["c1"]["cap_state"], [1.0]; atol=2e-1))
    end
    @testset "capcontrol_fbs" begin
        eng = parse_file("../test/data/opendss/IEEE13_CapControl.dss")
        sol = solve_mc_opf_capc(eng, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
        @test sol["termination_status"] == LOCALLY_SOLVED
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 0.404784; atol=5e-3)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]), -0.328146; atol=3e-1)
        @test all(isapprox.(sol["solution"]["bus"]["646"]["vm"], [1.03928, 1.05688]; atol=2e-1))
        @test all(isapprox.(sol["solution"]["bus"]["646"]["va"], [-89.997, 148.711]; atol=8e-1))
        @test all(isapprox.(sol["solution"]["shunt"]["c1"]["cap_state"], [1.0]; atol=6e-1))
    end
    @testset "capcontrol_lpubfdiag" begin
        eng = parse_file("../test/data/opendss/IEEE13_CapControl.dss")
        sol = solve_mc_opf_capc(eng, LPUBFDiagPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
        @test sol["termination_status"] == LOCALLY_SOLVED
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 0.404784; atol=5e-3)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]), -0.328146; atol=9e-1)
        @test all(isapprox.(sol["solution"]["bus"]["646"]["vm"], [1.03928, 1.05688]; atol=2e-1))
        @test all(isapprox.(sol["solution"]["shunt"]["c1"]["cap_state"], [1.0]; atol=6e-1))
    end
    @testset "capcontrol_fotr" begin
        eng = parse_file("../test/data/opendss/IEEE13_CapControl.dss")
        sol = solve_mc_opf_capc(eng, FOTRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
        @test sol["termination_status"] == LOCALLY_SOLVED
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 0.404784; atol=5e-3)
        @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]), -0.328146; atol=4e-1)
        @test all(isapprox.(sol["solution"]["bus"]["646"]["vm"], [1.03928, 1.05688]; atol=2e-1))
        @test all(isapprox.(sol["solution"]["bus"]["646"]["va"], [-89.997, 148.711]; atol=1e0))
        @test all(isapprox.(sol["solution"]["shunt"]["c1"]["cap_state"], [1.0]; atol=6e-1))
    end
end