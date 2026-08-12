@info "running multinetwork tests"

@testset "test multinetwork" begin
    @testset "3-bus balanced multinetwork nfa opb" begin
        math_ts = instantiate_mc_model_ravens(ravens_case3_balanced, NFAUPowerModel, build_mn_mc_opf; multinetwork=true)
        result_mn = optimize_model!(
            math_ts,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result_mn["termination_status"] == LOCALLY_SOLVED
    end

    @testset "3-bus balanced multinetwork instantiate_mc_model - autodetect multinetwork" begin
        pm_mn = instantiate_mc_model_ravens(ravens_case3_balanced, NFAUPowerModel, build_mn_mc_opf; multinetwork=true)
        @test ismultinetwork(pm_mn)
    end

    @testset "apply_voltage_bounds! to multinetworks" begin
        # math = transform_data_model()
        mn_math = instantiate_mc_model_ravens(ravens_case3_balanced, NFAUPowerModel, build_mn_mc_opf; multinetwork=true)
        
        apply_voltage_bounds_math!(mn_math)
        @debug keys(mn_math.data)
        for (n,nw) in mn_math.data["nw"]
            vbases, _ = calc_voltage_bases(mn_math.data["nw"]["1"], mn_math.data["nw"]["1"]["settings"]["vbases_default"])
            
            #TODO: figure out what is actually broken here
            @test all(all(isapprox.(bus["vmax"][filter(x->bus["grounded"][x] == 0,bus["terminals"])]/vbases[id], 1.1; atol=1e-6)) for (id,bus) in filter(x->x.second["name"]!="sourcebus",nw["bus"])) #TODO: vm == 1.1 but vbases != 1 thus the test fails
            @test all(all(isapprox.(bus["vmin"][filter(x->bus["grounded"][x] == 0,bus["terminals"])]/vbases[id], 0.9; atol=1e-6)) for (id,bus) in filter(x->x.second["name"]!="sourcebus",nw["bus"])) #TODO: vm ~= 0.9 but vbases != 1 thus the test fails  
        end
    end
    
    #IEEE13 Feeder does not convert from dss-->ravens
    # @testset "solve_mc_opf_oltc" begin
    #     math = transform_data_model(ravens_IEEE13_Feeder_engr)
    #     result_mn = PowerModelsDistribution.solve_mn_mc_opf_oltc(math, ACPUPowerModel, ipopt_solver)
    #     @test result_mn["termination_status"] == LOCALLY_SOLVED

    #     @test all(isapprox.(result_mn["solution"]["nw"]["1"]["voltage_source"]["source"]["pg"], [738.58786, 788.38272, 787.79729]; atol=1e-5))
    #     @test all(isapprox.(result_mn["solution"]["nw"]["1"]["voltage_source"]["source"]["qg"], [237.68517, 209.61208, 266.77223]; atol=1e-5))
    #     @test all(isapprox.(result_mn["solution"]["nw"]["8"]["voltage_source"]["source"]["pg"], [847.77707, 889.87745, 918.34146]; atol=1e-5))
    #     @test all(isapprox.(result_mn["solution"]["nw"]["8"]["voltage_source"]["source"]["qg"], [284.46267, 227.28860, 292.33564]; atol=1e-5))

    #     @test all(isapprox.(result_mn["solution"]["nw"]["1"]["transformer"]["reg1"]["tap"][2], [1.02358, 1.01724, 1.02169]; atol=1e-5))
    #     @test all(isapprox.(result_mn["solution"]["nw"]["8"]["transformer"]["reg1"]["tap"][2], [1.02719, 1.01984, 1.02414]; atol=1e-5))
    # end
end
