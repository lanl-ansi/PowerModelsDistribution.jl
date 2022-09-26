@info "running multinetwork tests"

@testset "test multinetwork" begin
    @testset "3-bus balanced multinetwork nfa opb" begin
        eng_ts = make_multinetwork(case3_balanced)
        result_mn = PowerModelsDistribution._solve_mn_mc_opb(eng_ts, NFAUPowerModel, ipopt_solver)

        @test result_mn["termination_status"] == LOCALLY_SOLVED
    end

    @testset "3-bus balanced multinetwork instantiate_mc_model - autodetect multinetwork" begin
        eng_ts = make_multinetwork(case3_balanced)
        pm_mn = instantiate_mc_model(eng_ts, ACPUPowerModel, build_mn_mc_opf)

        @test ismultinetwork(pm_mn)
    end

    @testset "apply_voltage_bounds! to multinetworks" begin
        mn_eng = make_multinetwork(case3_balanced)

        apply_voltage_bounds!(mn_eng)
        for (n,nw) in mn_eng["nw"]
            vbases, _ = calc_eng_voltage_bases(mn_eng["nw"]["1"], mn_eng["nw"]["1"]["settings"]["vbases_default"])

            @test all(all(isapprox.(bus["vm_ub"][filter(x->x∉bus["grounded"],bus["terminals"])]/vbases[id], 1.1; atol=1e-6)) for (id,bus) in filter(x->x.first!="sourcebus",nw["bus"]))
            @test all(all(isapprox.(bus["vm_lb"][filter(x->x∉bus["grounded"],bus["terminals"])]/vbases[id], 0.9; atol=1e-6)) for (id,bus) in filter(x->x.first!="sourcebus",nw["bus"]))
        end
    end
end
