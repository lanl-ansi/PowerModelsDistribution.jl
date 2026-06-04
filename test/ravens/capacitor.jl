@info "running capacitor control tests"

@testset "capacitor control" begin
    @testset "capcontrol_acp" begin
        math = transform_data_model(ravens_IEEE13_CapControl)
        result = solve_mc_opf_capc(math, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 405.556; atol=5)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), -527.15; atol=200)

        sourcebus_id = math["bus_lookup"]["646"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["646"]["vm"] ./ vbase, [1.03533, 1.06987]; atol=2e-1))
        @test all(isapprox.(result["bus"]["646"]["va"], [-90.1768, 148.498]; atol=3e-1))

        # @test all(isapprox.(result["shunt"]["c1"]["cap_state"], [1.0]; atol=2e-1))
    end

    @testset "capcontrol_acr" begin
        math = transform_data_model(ravens_IEEE13_CapControl)
        result = solve_mc_opf_capc(math, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == ALMOST_LOCALLY_SOLVED || result["termination_status"] == LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 405.556; atol=5)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), -527.15; atol=300)

        sourcebus_id = math["bus_lookup"]["646"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["646"]["vm"] ./ vbase, [1.03533, 1.06987]; atol=2e-1))
        @test all(isapprox.(result["bus"]["646"]["va"], [-90.1768, 148.498]; atol=3e-1))

        # @test all(isapprox.(result["shunt"]["c1"]["cap_state"], [1.0]; atol=2e-1))
    end

    @testset "capcontrol_ivr" begin
        math = transform_data_model(ravens_IEEE13_CapControl)
        result = solve_mc_opf_capc(math, IVRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 405.556; atol=5)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), -527.15; atol=300)

        sourcebus_id = math["bus_lookup"]["646"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["646"]["vm"] ./ vbase, [1.03533, 1.06987]; atol=2e-1))
        @test all(isapprox.(result["bus"]["646"]["va"], [-90.1768, 148.498]; atol=3e-1))

        # @test all(isapprox.(result["shunt"]["c1"]["cap_state"], [1.0]; atol=2e-1))
    end

    @testset "capcontrol_fbs" begin
        math = transform_data_model(ravens_IEEE13_CapControl)
        result = solve_mc_opf_capc(IEEE13_CapControl, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["pg"]), 404.784; atol=5)
        @test isapprox(sum(result["solution"]["voltage_source"]["source"]["qg"]), -328.146; atol=300)

        sourcebus_id = math["bus_lookup"]["646"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["solution"]["bus"]["646"]["vm"] ./ vbase, [1.03928, 1.05688]; atol=2e-1))
        @test all(isapprox.(result["solution"]["bus"]["646"]["va"], [-89.997, 148.711]; atol=8e-1))

        # @test all(isapprox.(result["shunt"]["c1"]["cap_state"], [1.0]; atol=6e-1))
    end

    @testset "capcontrol_lpubfdiag" begin
        math = transform_data_model(ravens_IEEE13_CapControl)
        result = solve_mc_opf_capc(math, LPUBFDiagPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 404.784; atol=5)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), -328.146; atol=900)

        sourcebus_id = math["bus_lookup"]["646"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["646"]["vm"] ./ vbase, [1.03928, 1.05688]; atol=2e-1))

        # @test all(isapprox.(result["shunt"]["c1"]["cap_state"], [1.0]; atol=6e-1))
    end

    @testset "capcontrol_fotr" begin
        math = transform_data_model(ravens_IEEE13_CapControl)
        result = solve_mc_opf_capc(math, FOTRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

        result = transform_solution(result["solution"], math, make_si=true)

        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 404.784; atol=5)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), -328.146; atol=400)
    
        sourcebus_id = math["bus_lookup"]["646"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["646"]["vm"] ./ vbase, [1.03928, 1.05688]; atol=2e-1))
        @test all(isapprox.(result["bus"]["646"]["va"], [-89.997, 148.711]; atol=1e0))

        # @test all(isapprox.(result["shunt"]["c1"]["cap_state"], [1.0]; atol=6e-1))
    end

    @testset "capcontrol_fotp" begin
        math = transform_data_model(ravens_IEEE13_CapControl)
        result = solve_mc_opf_capc(math, FOTPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED


        result = transform_solution(result["solution"], math, make_si=true)


        @test isapprox(sum(result["voltage_source"]["source"]["pg"]), 404.784; atol=5)
        @test isapprox(sum(result["voltage_source"]["source"]["qg"]), -328.146; atol=400)

        sourcebus_id = math["bus_lookup"]["646"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]

        @test all(isapprox.(result["bus"]["646"]["vm"] ./ vbase, [1.03928, 1.05688]; atol=2e-1))
        @test all(isapprox.(result["bus"]["646"]["va"], [-89.997, 148.711]; atol=1e0))

        # @test all(isapprox.(result["shunt"]["c1"]["cap_state"], [1.0]; atol=6e-1))
    end
end
