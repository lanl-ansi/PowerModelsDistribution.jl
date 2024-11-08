@info "running optimal power flow (opf) tests using RAVENS data"

@testset "test opf ravens" begin

    @testset "ravens case 3 with gens" begin
        pmd_model = instantiate_mc_model_ravens(ravens_case3_withgens, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED
    end

    @testset "ravens case 3 with PV and storage" begin
        pmd_model = instantiate_mc_model_ravens(ravens_case3_withpvandstorage, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED
    end

    @testset "ravens case 3 with subxf" begin
        pmd_model = instantiate_mc_model_ravens(ravens_case3_withsubxf, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED
    end

    @testset "ravens case 3 with capacitor" begin
        pmd_model = instantiate_mc_model_ravens(ravens_case3_withcap, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED
    end

    @testset "ravens test with switches 3w" begin
        pmd_model = instantiate_mc_model_ravens(ravens_test_switch_3w, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED
    end

    # @testset "ravens test with switches 1w" begin
    #     pmd_model = instantiate_mc_model_ravens(ravens_test_switch_1w, ACPUPowerModel, build_mc_opf)
    #     result = optimize_model!(
    #         pmd_model,
    #         relax_integrality=false,
    #         optimizer=ipopt_solver,
    #         solution_processors=Function[]
    #     )
    #     @test result["termination_status"] == LOCALLY_SOLVED
    # end
end
