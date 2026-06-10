@info "running power flow (pf) tests - extended to RAVENS"

#TODO: do we want all of the conversions [/4.33, *10^3, radians to degrees]

@testset "test pf - ravens" begin
    @testset "2-bus diagonal acp pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case2_diag), ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        math = transform_data_model(ravens_case2_diag)
        primary_id = math["bus_lookup"]["primary"]

        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["vm"]./4.33012636656, 0.227339; atol=1e-4))

        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["va"], deg2rad.([-0.657496, -120.657, 119.343]); atol=1e-2))

        @test isapprox(sum(result["solution"]["gen"]["1"]["pg"])*10^3, 18.20887; atol=1e-4)
        @test isapprox(sum(result["solution"]["gen"]["1"]["qg"])*10^3,  0.20887; atol=1e-4)
    end

    @testset "2-bus diagonal acr pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case2_diag), ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        math = transform_data_model(ravens_case2_diag)
        primary_id = math["bus_lookup"]["primary"]

        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["vm"]./4.33012636656, 0.227339; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["va"], deg2rad.([-0.657496, -120.657, 119.343]); atol=1e-2))

        @test isapprox(sum(result["solution"]["gen"]["1"]["pg"])*10^3, 18.20887; atol=1e-4)
        @test isapprox(sum(result["solution"]["gen"]["1"]["qg"])*10^3,  0.20887; atol=1e-4)
    end

    @testset "2-bus diagonal ivr pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case2_diag), IVRUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        @test isapprox(sum(result["solution"]["gen"]["1"]["pg"])*10^3, 18.20896; atol=1e-4)
        @test isapprox(sum(result["solution"]["gen"]["1"]["qg"])*10^3,  0.20896; atol=1e-4)
    end

    @testset "3-bus balanced acp pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case3_balanced), ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        math = transform_data_model(ravens_case3_balanced)
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        primary_id = math["bus_lookup"]["primary"]
        loadbus_id = math["bus_lookup"]["loadbus"]


        @test all(isapprox.(result["solution"]["bus"][string(sourcebus_id)]["vm"]./4.33012636656, 0.229993; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(sourcebus_id)]["va"], deg2rad.([0.0, -120.0, 120.0]); atol=1e-2))

        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["vm"]./4.33012636656, 0.227932; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["va"], deg2rad.([-0.08, -120.08, 119.92]); atol=0.2))

        @test all(isapprox.(result["solution"]["bus"][string(loadbus_id)]["vm"]./4.33012636656, 0.225537; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(loadbus_id)]["va"], deg2rad.([-0.17, -120.17, 119.83]); atol=0.2))

        @test isapprox(sum(result["solution"]["gen"]["1"]["pg"])*10^3, 18.34478; atol=1e-2)
        @test isapprox(sum(result["solution"]["gen"]["1"]["qg"])*10^3,  9.19392; atol=1e-2)
    end

    @testset "3-bus balanced acr pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case3_balanced), ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        math = transform_data_model(ravens_case3_balanced)
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        primary_id = math["bus_lookup"]["primary"]
        loadbus_id = math["bus_lookup"]["loadbus"]

        @test all(isapprox.(result["solution"]["bus"][string(sourcebus_id)]["vm"]./4.33012636656, 0.229993; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(sourcebus_id)]["va"], deg2rad.([0.0, -120.0, 120.0]); atol=1e-2))

        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["vm"]./4.33012636656, 0.227932; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["va"], deg2rad.([-0.08, -120.08, 119.92]); atol=0.2))

        @test all(isapprox.(result["solution"]["bus"][string(loadbus_id)]["vm"]./4.33012636656, 0.225537; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(loadbus_id)]["va"], deg2rad.([-0.17, -120.17, 119.83]); atol=0.2))

        @test isapprox(sum(result["solution"]["gen"]["1"]["pg"])*10^3, 18.34478; atol=1e-2)
        @test isapprox(sum(result["solution"]["gen"]["1"]["qg"])*10^3,  9.19392; atol=1e-2)
    end

    @testset "3-bus balanced ivr pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case3_balanced), IVRUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        @test isapprox(sum(result["solution"]["gen"]["1"]["pg"])*10^3, 18.34498; atol=1e-5)
        @test isapprox(sum(result["solution"]["gen"]["1"]["qg"])*10^3,  9.19404; atol=1e-4)
    end

    @testset "3-bus balanced no linecode basefreq defined acp pf" begin
        result1 = solve_mc_pf(transform_data_model(ravens_case3_balanced), ACPUPowerModel, ipopt_solver)
        result2 = solve_mc_pf(transform_data_model(ravens_case3_balanced_basefreq), ACPUPowerModel, ipopt_solver)

        math1 = transform_data_model(ravens_case3_balanced)
        sourcebus_id1 = math1["bus_lookup"]["sourcebus"]
        primary_id1 = math1["bus_lookup"]["primary"]
        loadbus_id1 = math1["bus_lookup"]["loadbus"]

        math2 = transform_data_model(ravens_case3_balanced_basefreq)
        sourcebus_id2 = math2["bus_lookup"]["sourcebus"]
        primary_id2 = math2["bus_lookup"]["primary"]
        loadbus_id2 = math2["bus_lookup"]["loadbus"]

        @test all(all(isapprox.(result1["solution"]["bus"][string(sourcebus_id1)]["vm"]./4.33012636656, result2["solution"]["bus"][string(sourcebus_id1)]["vm"]./4.33012636656; atol=1e-4)) for (i, bus) in result1["solution"]["bus"])
        @test all(all(isapprox.(result1["solution"]["bus"][string(sourcebus_id1)]["va"], result2["solution"]["bus"][string(sourcebus_id1)]["va"]; atol=1e-4)) for (i, bus) in result1["solution"]["bus"])

        @test all(all(isapprox.(result1["solution"]["bus"][string(primary_id1)]["vm"]./4.33012636656, result2["solution"]["bus"][string(primary_id1)]["vm"]./4.33012636656; atol=1e-4)) for (i, bus) in result1["solution"]["bus"])
        @test all(all(isapprox.(result1["solution"]["bus"][string(primary_id1)]["va"], result2["solution"]["bus"][string(primary_id1)]["va"]; atol=1e-4)) for (i, bus) in result1["solution"]["bus"])

        @test all(all(isapprox.(result1["solution"]["bus"][string(loadbus_id1)]["vm"]./4.33012636656, result2["solution"]["bus"][string(loadbus_id1)]["vm"]./4.33012636656; atol=1e-4)) for (i, bus) in result1["solution"]["bus"])
        @test all(all(isapprox.(result1["solution"]["bus"][string(loadbus_id1)]["va"], result2["solution"]["bus"][string(loadbus_id1)]["va"]; atol=1e-4)) for (i, bus) in result1["solution"]["bus"])

        @test isapprox(sum(result1["solution"]["gen"]["1"]["pg"])*10^3, sum(result2["solution"]["gen"]["1"]["pg"])*10^3; atol=1e-4)
        @test isapprox(sum(result1["solution"]["gen"]["1"]["qg"])*10^3, sum(result2["solution"]["gen"]["1"]["qg"])*10^3; atol=1e-2)
    end

    @testset "3-bus unbalanced acp pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case3_unbalanced), ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
       
        math = transform_data_model(ravens_case3_unbalanced)
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        vbase = vbase * 4.33012636656

        primary_id = math["bus_lookup"]["primary"]
        loadbus_id = math["bus_lookup"]["loadbus"]

        @test all(isapprox.(result["solution"]["bus"][string(sourcebus_id)]["vm"] ./ vbase, [0.9959, 0.9959, 0.9959]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(sourcebus_id)]["va"], deg2rad.([0.0, -120.0, 120.0]); atol=1e-2))

        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["vm"] ./ vbase, [0.98094, 0.989365, 0.987043]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["va"], deg2rad.([-0.22, -120.11, 120.12]); atol=1e-2))

        @test all(isapprox.(result["solution"]["bus"][string(loadbus_id)]["vm"] ./ vbase, [0.96355, 0.981767, 0.976786]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(loadbus_id)]["va"], deg2rad.([-0.48, -120.24, 120.27]); atol=1e-2))

        @test isapprox(sum(result["solution"]["gen"]["1"]["pg"])*10^3, 21.4812; atol=1e-2)
        @test isapprox(sum(result["solution"]["gen"]["1"]["qg"])*10^3, 9.27263; atol=1e-2)
    end

    @testset "3-bus unbalanced acr pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case3_unbalanced), ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        math = transform_data_model(ravens_case3_unbalanced)
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        vbase = vbase * 4.33012636656

        primary_id = math["bus_lookup"]["primary"]
        loadbus_id = math["bus_lookup"]["loadbus"]

        @test all(isapprox.(result["solution"]["bus"][string(sourcebus_id)]["vm"] ./ vbase, [0.9959, 0.9959, 0.9959]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(sourcebus_id)]["va"], deg2rad.([0.0, -120.0, 120.0]); atol=1e-2))

        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["vm"] ./ vbase, [0.98094, 0.989365, 0.987043]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["va"], deg2rad.([-0.22, -120.11, 120.12]); atol=1e-2))

        @test all(isapprox.(result["solution"]["bus"][string(loadbus_id)]["vm"] ./ vbase, [0.96355, 0.981767, 0.976786]; atol=1e-4))
        @test all(isapprox.(result["solution"]["bus"][string(loadbus_id)]["va"], deg2rad.([-0.48, -120.24, 120.27]); atol=1e-2))

        @test isapprox(sum(result["solution"]["gen"]["1"]["pg"])*10^3, 21.4812; atol=1e-2)
        @test isapprox(sum(result["solution"]["gen"]["1"]["qg"])*10^3, 9.27263; atol=1e-2)
    end

    @testset "3-bus unbalanced w/ asymmetric linecode & phase order swap acp pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case3_unbalanced_assym_swap), ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        math = transform_data_model(ravens_case3_unbalanced_assym_swap)
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        vbase = vbase * 4.33012636656

        primary_id = math["bus_lookup"]["primary"]
        
        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["vm"] ./ vbase, [0.983453, 0.98718, 0.981602]; atol=1e-4))

        @test all(isapprox.(result["solution"]["bus"][string(primary_id)]["va"], deg2rad.([-0.07, -120.19, 120.29]); atol=1e-2))
    end

    @testset "5-bus phase drop acp pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case5_phase_drop), ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        math = transform_data_model(ravens_case5_phase_drop)
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        vbase = vbase * 4.33012636656

        
        midbus_id = math["bus_lookup"]["midbus"]
        @test all(isapprox.(result["solution"]["bus"][string(midbus_id)]["vm"] ./ vbase, [0.973519, 0.964902, 0.956465]; atol = 1e-4))
    end

    @testset "5-bus phase drop acr pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case5_phase_drop), ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        @test result["termination_status"] == LOCALLY_SOLVED

        math = transform_data_model(ravens_case5_phase_drop)
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        vbase = vbase * 4.33012636656

        
        midbus_id = math["bus_lookup"]["midbus"]
        @test all(isapprox.(result["solution"]["bus"][string(midbus_id)]["vm"] ./ vbase, [0.973519, 0.964902, 0.956465]; atol=1e-4))
    end

    @testset "matrix branch shunts acp pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case2_mxshunt), ACPUPowerModel, ipopt_solver)

        math = transform_data_model(ravens_case2_mxshunt)
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        vbase = vbase * 4.33012636656


        
        loadbus_id = math["bus_lookup"]["loadbus"]
        @test all(isapprox.(result["solution"]["bus"][string(loadbus_id)]["vm"] ./ vbase, [0.987399, 0.981300, 1.003536]; atol=1e-4))
    end

    @testset "matrix branch shunts acr pf" begin
        result = solve_mc_pf(transform_data_model(ravens_case2_mxshunt), ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        math = transform_data_model(ravens_case2_mxshunt)
        sourcebus_id = math["bus_lookup"]["sourcebus"]
        vbase = math["bus"][string(sourcebus_id)]["vbase"]
        vbase = vbase * 4.33012636656


        loadbus_id = math["bus_lookup"]["loadbus"]
        @test all(isapprox.(result["solution"]["bus"][string(loadbus_id)]["vm"] ./ vbase, [0.987399, 0.981300, 1.003536]; atol=1e-4))
    end

    @testset "virtual sourcebus creation acp pf" begin
        result = solve_mc_pf(transform_data_model(case2_virtual_sourcebus), ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED

        @test all(all(isapprox.(result["solution"]["bus"]["$n"]["vm"], [0.961352, 0.999418, 1.00113]; atol=1e-4)) for n in [1, 2])
        @test all(all(isapprox.(result["solution"]["bus"]["$n"]["va"], deg2rad.([-1.25, -120.06, 120.0]); atol=1e-2)) for n in [1, 2])
    end
end
