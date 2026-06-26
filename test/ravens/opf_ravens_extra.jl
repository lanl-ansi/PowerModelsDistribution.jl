@info "running optimal power flow (opf) tests using extra RAVENS data"

@testset "test opf ravens extra" begin

    @testset "Test case2 virtual sourcebus" begin
        pmd_model = instantiate_mc_model_ravens(ravens_case2_virtual_sourcebus, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
    end

    @testset "Test case3 balanced battery" begin
        pmd_model = instantiate_mc_model_ravens(ravens_case3_balanced_battery, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
    end

    @testset "Test case3 lm models" begin
        pmd_model = instantiate_mc_model_ravens(ravens_case3_lm_models, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_INFEASIBLE
    end

    @testset "Test trans 2w yy" begin
        pmd_model = instantiate_mc_model_ravens(ravens_ut_trans_2w_yy, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Test trans 2w dy" begin
        pmd_model = instantiate_mc_model_ravens(ravens_ut_trans_2w_dy_lead, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Test trans 2w yy lag" begin
        pmd_model = instantiate_mc_model_ravens(ravens_ut_trans_2w_dy_lag, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Test trans 2w yy lssi" begin
        pmd_model = instantiate_mc_model_ravens(ravens_ut_trans_2w_dy_lead_small_series_impedance, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED
    end

    @testset "Test trans 2w yy bank" begin
        pmd_model = instantiate_mc_model_ravens(ravens_ut_trans_2w_yy_bank, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_INFEASIBLE
    end

    @testset "Test trans 2w yy unbanked" begin
        pmd_model = instantiate_mc_model_ravens(ravens_ut_trans_2w_yy_unbanked, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_INFEASIBLE
    end

    @testset "Test trans 3w dyy 1" begin
        pmd_model = instantiate_mc_model_ravens(ravens_ut_trans_3w_dyy_1, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_INFEASIBLE
    end

    @testset "Test trans 3w dyy 2" begin
        pmd_model = instantiate_mc_model_ravens(ravens_ut_trans_3w_dyy_2, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_INFEASIBLE
    end

    @testset "Test trans 3w dyy 3" begin
        pmd_model = instantiate_mc_model_ravens(ravens_ut_trans_3w_dyy_3, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_INFEASIBLE
    end

    @testset "Test trans 3w dyy 3 loadloss" begin
        pmd_model = instantiate_mc_model_ravens(ravens_ut_trans_3w_dyy_3_loadloss, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_INFEASIBLE
    end

    @testset "Test trans 3w center tap" begin
        @test_broken begin
            pmd_model = instantiate_mc_model_ravens(ravens_trans_3w_center_tap, ACPUPowerModel, build_mc_opf)
            result = optimize_model!(
                pmd_model,
                relax_integrality=false,
                optimizer=ipopt_solver,
                solution_processors=Function[]
            )
            @test result["termination_status"] == LOCALLY_INFEASIBLE
        end
    end

    @testset "Test IEEE13 Assets" begin
        pmd_model = instantiate_mc_model_ravens(ravens_IEEE13_Assets, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED #LOCALLY_INFEASIBLE
    end

    @testset "Test IEEE13 Reg Control" begin
        pmd_model = instantiate_mc_model_ravens(ravens_IEEE13_RegControl, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
    end


    @testset "Test IEEE13 Cap Control" begin
        pmd_model = instantiate_mc_model_ravens(ravens_IEEE13_CapControl, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
    end

    @testset "Test gen 3ph wye" begin
        pmd_model = instantiate_mc_model_ravens(ravens_test_gen_3ph_wye, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
    end

    @testset "Test switch" begin
        pmd_model = instantiate_mc_model_ravens(ravens_test_switch, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED
    end

    @testset "Test trans dy" begin
        pmd_model = instantiate_mc_model_ravens(ravens_test_trans_dy, ACPUPowerModel, build_mc_opf)
        result = optimize_model!(
            pmd_model,
            relax_integrality=false,
            optimizer=ipopt_solver,
            solution_processors=Function[]
        )
        @test result["termination_status"] == LOCALLY_SOLVED
    end
    

end
