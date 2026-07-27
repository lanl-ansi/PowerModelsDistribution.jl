@info "running switch tests"

eng = PowerModelsDistribution.parse_file("../test/data/opendss/case3_balanced_switch.dss")
eng["switch"]["ohline"]["length"] = 1.0

@testset "test switch opf" begin
    @testset "3-bus balanced switch closed opf acp" begin
       result = solve_mc_opf(eng, ACPUPowerModel, ipopt_solver_adaptive; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"], [0.0, -0.032919, -0.0663864], [0.9959, 0.986973, 0.976606])
            @test all(isapprox.(result["solution"]["bus"][bus]["va"], [0, -120, 120] .+ va; atol=0.2))
            @test all(isapprox.(result["solution"]["bus"][bus]["vm"], vm; atol=1e-3))
        end
    end

    @testset "3-bus balanced switch closed opf acr" begin
        result = solve_mc_opf(eng, ACRUPowerModel, ipopt_solver_adaptive; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"], [0.0, -0.032919, -0.066387], [0.995899, 0.9869727, 0.9766053])
            @test all(isapprox.(calc_va_acr(result, bus), [0, -120, 120] .+ va; atol=0.2))
            @test all(isapprox.(calc_vm_acr(result, bus), vm; atol=1e-3))
        end
    end

    @testset "3-bus balanced switch closed opf ivr" begin
        result = solve_mc_opf(eng, IVRUPowerModel, ipopt_solver_adaptive; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"], [0.0, -0.032919, -0.066387], [0.995899, 0.9869727, 0.9766053])
            @test all(isapprox.(calc_va_acr(result, bus), [0, -120, 120] .+ va; atol=0.2))
            @test all(isapprox.(calc_vm_acr(result, bus), vm; atol=1e-3))
        end
    end

    @testset "3-bus balanced switch closed opf nfa" begin
        result = solve_mc_opf(eng, NFAUPowerModel, ipopt_solver_adaptive; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test all(isapprox.(result["solution"]["voltage_source"]["source"]["pg"], 6.0e-5; atol=1e-5))
    end

    @testset "3-bus balanced switch closed opf lpubfdiag" begin
        result = solve_mc_opf(eng, LPUBFDiagPowerModel, ipopt_solver_adaptive; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        for (bus, vm) in zip(["sourcebus", "primary", "loadbus"], [0.9959, 0.987107, 0.976797])
            @test all(isapprox.(calc_vm_w(result, bus), vm; atol=1e-3))
        end
    end
end

@testset "test switch pf" begin
    @testset "3-bus balanced switch closed pf acp" begin
        result = solve_mc_pf(eng, ACPUPowerModel, ipopt_solver_adaptive; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED

        for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"], [0.0, -0.032919, -0.0663864], [0.9959, 0.986973, 0.976606])
            @test all(isapprox.(result["solution"]["bus"][bus]["va"], [0, -120, 120] .+ va; atol=0.2))
            @test all(isapprox.(result["solution"]["bus"][bus]["vm"], vm; atol=1e-3))
        end
    end

    @testset "3-bus balanced switch closed pf acr" begin
        result = solve_mc_pf(eng, ACRUPowerModel, ipopt_solver_adaptive; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED

        for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"], [0.0, -0.032919, -0.066387], [0.995899, 0.9869727, 0.9766053])
            @test all(isapprox.(calc_va_acr(result, bus), [0, -120, 120] .+ va; atol=0.2))
            @test all(isapprox.(calc_vm_acr(result, bus), vm; atol=1e-3))
        end
    end

    @testset "3-bus balanced switch closed pf ivr" begin
        result = solve_mc_pf(eng, IVRUPowerModel, ipopt_solver_adaptive; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"], [0.0, -0.032919, -0.066387], [0.995899, 0.9869727, 0.9766053])
            @test all(isapprox.(calc_va_acr(result, bus), [0, -120, 120] .+ va; atol=0.2))
            @test all(isapprox.(calc_vm_acr(result, bus), vm; atol=1e-3))
        end
    end

    @testset "3-bus balanced switch closed pf nfa" begin
        result = solve_mc_pf(eng, NFAUPowerModel, ipopt_solver_adaptive; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test all(isapprox.(result["solution"]["voltage_source"]["source"]["pg"], 6.0e-5; atol=1e-6))
    end

    @testset "3-bus balanced switch closed pf lpubfdiag" begin
        result = solve_mc_pf(eng, LPUBFDiagPowerModel, ipopt_solver_adaptive; make_si=false)

        @test result["termination_status"] == LOCALLY_SOLVED
        for (bus, vm) in zip(["sourcebus", "primary", "loadbus"], [0.9959, 0.987107, 0.976797])
            @test all(isapprox.(calc_vm_w(result, bus), vm; atol=1e-3))
        end
    end

end

@testset "test switch mld" begin
    eng["switch"]["ohline"]["state"] = OPEN

    @testset "3-bus balanced switch open mld acp" begin
        result = solve_mc_mld(eng, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        for (_,load) in result["solution"]["load"]
            @test all(isapprox.(load["status"], 0.0; atol=1e-5))
        end
    end

    @testset "3-bus balanced switch open mld acr" begin
        result = solve_mc_mld(eng, ACRUPowerModel, ipopt_solver_adaptive)

        @test result["termination_status"] == LOCALLY_SOLVED
        for (_,load) in result["solution"]["load"]
            @test all(isapprox.(load["status"], 0.0; atol=1e-5))
        end
    end

    @testset "3-bus balanced switch open mld nfa" begin
        result = solve_mc_mld(eng, NFAUPowerModel, ipopt_solver_adaptive)

        @test result["termination_status"] == LOCALLY_SOLVED
        for (_,load) in result["solution"]["load"]
            @test all(isapprox.(load["status"], 0.0; atol=1e-5))
        end
    end

    @testset "3-bus balanced switch open mld lpubfdiag" begin
        result = solve_mc_mld(eng, LPUBFDiagPowerModel, ipopt_solver_adaptive)

        @test result["termination_status"] == LOCALLY_SOLVED
        for (_,load) in result["solution"]["load"]
            @test all(isapprox.(load["status"], 0.0; atol=1e-5))
        end
    end
end

@testset "test switch make lossless" begin
    eng["switch"]["ohline"]["state"] = CLOSED
    make_lossless!(eng)

    result = solve_mc_opf(eng, ACPUPowerModel, ipopt_solver_adaptive)

    @test result["termination_status"] == LOCALLY_SOLVED
    @test all(result["solution"]["switch"]["ohline"]["psw_fr"] .== -result["solution"]["switch"]["ohline"]["psw_to"])
    @test all(result["solution"]["switch"]["ohline"]["qsw_fr"] .== -result["solution"]["switch"]["ohline"]["qsw_to"])
end
