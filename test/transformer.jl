@info "running transformer tests"

@testset "transformers" begin
    @testset "test transformer acp pf" begin
        @testset "2w transformer acp pf yy" begin
            eng = parse_file("../test/data/opendss/ut_trans_2w_yy.dss")
            sol = solve_mc_pf(eng, ACPUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.87451, 0.8613, 0.85348], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[-0.1, -120.4, 119.8], Inf) <= 0.1
        end

        @testset "2w transformer acp pf dy_lead" begin
            eng = parse_file("../test/data/opendss/ut_trans_2w_dy_lead.dss")
            sol = solve_mc_pf(eng, ACPUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.87391, 0.86055, 0.85486], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[29.8, -90.4, 149.8], Inf) <= 0.1
        end

        @testset "2w transformer acp pf dy_lag" begin
            eng = parse_file("../test/data/opendss/ut_trans_2w_dy_lag.dss")
            sol = solve_mc_pf(eng, ACPUPowerModel, ipopt_solver; make_si=false)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.92092, 0.91012, 0.90059], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[-30.0, -150.4, 89.8], Inf) <= 0.1
        end
    end

    @testset "test transformer ivr pf" begin
        @testset "2w transformer ivr pf yy" begin
            eng = parse_file("../test/data/opendss/ut_trans_2w_yy.dss")
            sol = solve_mc_pf(eng, IVRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.87451, 0.8613, 0.85348], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[-0.1, -120.4, 119.8], Inf) <= 0.1
        end

        @testset "2w transformer ivr pf dy_lead" begin
            eng = parse_file("../test/data/opendss/ut_trans_2w_dy_lead.dss")
            sol = solve_mc_pf(eng, IVRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.87391, 0.86055, 0.85486], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[29.8, -90.4, 149.8], Inf) <= 0.1
        end

        @testset "2w transformer ivr pf dy_lag" begin
            eng = parse_file("../test/data/opendss/ut_trans_2w_dy_lag.dss")
            sol = solve_mc_pf(eng, IVRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.92092, 0.91012, 0.90059], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[-30.0, -150.4, 89.8], Inf) <= 0.1
        end
    end

    @testset "test transformer acr pf" begin
        @testset "2w transformer acr pf yy" begin
            eng = parse_file("../test/data/opendss/ut_trans_2w_yy.dss")
            sol = solve_mc_pf(eng, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.87451, 0.8613, 0.85348], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[-0.1, -120.4, 119.8], Inf) <= 0.1
        end

        @testset "2w transformer acr pf dy_lead" begin
            eng = parse_file("../test/data/opendss/ut_trans_2w_dy_lead.dss")
            sol = solve_mc_pf(eng, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.87391, 0.86055, 0.85486], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[29.8, -90.4, 149.8], Inf) <= 0.1
        end

        @testset "2w transformer acr pf dy_lag" begin
            eng = parse_file("../test/data/opendss/ut_trans_2w_dy_lag.dss")
            sol = solve_mc_pf(eng, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.92092, 0.91012, 0.90059], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[-30.0, -150.4, 89.8], Inf) <= 0.1
        end
    end

    @testset "2w transformer ac pf yy - banked transformers" begin
        eng1 = parse_file("../test/data/opendss/ut_trans_2w_yy_bank.dss")
        eng2 = parse_file("../test/data/opendss/ut_trans_2w_yy_bank.dss"; bank_transformers=false)
        result1 = solve_mc_pf(eng1, ACPUPowerModel, ipopt_solver)
        result2 = solve_mc_pf(eng2, ACPUPowerModel, ipopt_solver)

        @test result1["termination_status"] == LOCALLY_SOLVED
        @test result2["termination_status"] == LOCALLY_SOLVED

        # @test result1["solution"]["bus"] == result2["solution"]["bus"] # TODO need a new test, transformer model changed, use voltages on real bus
        # @test result1["solution"]["gen"] == result2["solution"]["gen"] # TODO need a new test, transformer model changed, use voltages on real bus
    end

    @testset "three winding transformer pf" begin
        @testset "3w transformer ac pf dyy - all non-zero"  begin
            file =
            eng = parse_file("../test/data/opendss/ut_trans_3w_dyy_1.dss")
            sol = solve_mc_pf(eng, ACPUPowerModel, ipopt_solver; make_si=false)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.9318, 0.88828, 0.88581], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[30.1, -90.7, 151.2], Inf) <= 0.1
        end

        @testset "3w transformer ac pf dyy - some non-zero" begin
            file =
            eng = parse_file("../test/data/opendss/ut_trans_3w_dyy_2.dss")
            sol = solve_mc_pf(eng, ACPUPowerModel, ipopt_solver; make_si=false)
            #@test isapprox(vm(sol, eng, "3"), [0.93876, 0.90227, 0.90454], atol=1E-4)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.93876, 0.90227, 0.90454], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[31.6, -88.8, 153.3], Inf) <= 0.1
        end

        @testset "3w transformer ac pf dyy - all zero" begin
            eng = parse_file("../test/data/opendss/ut_trans_3w_dyy_3.dss")
            sol = solve_mc_pf(eng, ACPUPowerModel, ipopt_solver; make_si=false)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.97047, 0.93949, 0.946], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[30.6, -90.0, 151.9], Inf) <= 0.1
        end

        @testset "3w transformer ac pf dyy - %loadloss=0" begin
            eng = parse_file("../test/data/opendss/ut_trans_3w_dyy_3_loadloss.dss")
            sol = solve_mc_pf(eng, ACPUPowerModel, ipopt_solver; make_si=false)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.969531, 0.938369, 0.944748], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-[30.7, -90.0, 152.0], Inf) <= 0.1
        end
    end

    @testset "oltc tests" begin
        @testset "2w transformer acp opf_oltc yy" begin
            eng = parse_file("../test/data/opendss/ut_trans_2w_yy_oltc.dss")

            # free the taps
            eng["transformer"]["tx1"]["tm_fix"] = fill(zeros(Bool, 3), 2)

            math = transform_data_model(eng)
            pm = instantiate_mc_model(math, ACPUPowerModel, build_mc_opf_oltc)
            sol = PMD.optimize_model!(pm, optimizer=ipopt_solver)

            # check that taps are set as to boost the voltage in the branches as much as possible;
            # this is trivially optimal if the voltage bounds are not binding
            # and without significant shunts (both branch and transformer)
            @test norm(tap(1,pm)-[0.95, 0.95, 0.95], Inf) <= 1E-4
            @test norm(tap(2,pm)-[1.05, 1.05, 1.05], Inf) <= 1E-4

            # then check whether voltage is what OpenDSS expects for those values
            solution = transform_solution(sol["solution"], math, make_si=false)

            @test norm(solution["bus"]["3"]["vm"]-[1.0352, 1.022, 1.0142], Inf) <= 1E-4
            @test norm(solution["bus"]["3"]["va"]-[-0.1, -120.4, 119.8], Inf) <= 0.1
        end
    end
    @testset "linearized transformers" begin
        @testset "2w_dy_lead" begin
            eng = parse_file("../test/data/opendss/ut_trans_2w_dy_lead.dss")
            sol = solve_mc_opf(eng, LPUBFDiagPowerModel, ipopt_solver)
            @test norm(sol["solution"]["bus"]["3"]["w"]-[0.76674, 0.74840, 0.73846], Inf) <= 1E-4
        end
        @testset "3w_dyy_1" begin
            eng = parse_file("../test/data/opendss/ut_trans_3w_dyy_1.dss")
            sol = solve_mc_opf(eng, LPUBFDiagPowerModel, ipopt_solver)
            @test norm(sol["solution"]["bus"]["3"]["w"]-[0.86095, 0.81344, 0.80480], Inf) <= 1E-4
        end
        @testset "3w_dyy_2" begin
            eng = parse_file("../test/data/opendss/ut_trans_3w_dyy_2.dss")
            sol = solve_mc_opf(eng, LPUBFDiagPowerModel, ipopt_solver)
            @test norm(sol["solution"]["bus"]["3"]["w"]-[0.87086, 0.83270, 0.83208], Inf) <= 1E-4
        end
        @testset "2w_dy_lead_small_series_impedance" begin
            eng = parse_file("../test/data/opendss/ut_trans_2w_dy_lead_small_series_impedance.dss", data_model=MATHEMATICAL)
            sola = solve_mc_opf(eng, LPUBFDiagPowerModel, ipopt_solver)
            solb = solve_mc_opf(eng, ACPUPowerModel, ipopt_solver)
            @test norm(sola["solution"]["bus"]["1"]["w"]-solb["solution"]["bus"]["1"]["vm"].^2, Inf) <= 1.2E-3
            @test norm(sola["solution"]["branch"]["1"]["pf"]-solb["solution"]["branch"]["1"]["pf"], Inf) <= 1E-3
        end
    end
    @testset "voltage regulator control" begin
        @testset "regcontrol_acp" begin
            eng = parse_file("../test/data/opendss/IEEE13_RegControl.dss")
            math = transform_data_model(eng)
            sol = solve_mc_opf_oltc(math, ACPUPowerModel, ipopt_solver; make_si=false)
            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"]), 0.405235; atol=5e-4)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"]), 0.2666; atol=9e-4)
            @test all(isapprox.(sol["solution"]["bus"]["1"]["vm"], [1.0189, 1.0313, 1.0313]; atol=2e-2))
            @test all(isapprox.(sol["solution"]["bus"]["1"]["va"]*180/pi, [30, -90, 150]; atol=3e-3))
            @test all(isapprox.(sol["solution"]["transformer"]["2"]["tap"], [1.01875, 1.03125, 1.03125]; atol=2e-2))
        end
        @testset "regcontrol_acr" begin
            eng = parse_file("../test/data/opendss/IEEE13_RegControl.dss")
            math = transform_data_model(eng)
            sol = solve_mc_opf_oltc(math, ACRUPowerModel, ipopt_solver; make_si=false, solution_processors=[sol_data_model!])
            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"]), 0.405235; atol=5e-4)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"]), 0.2666; atol=9e-4)
            @test all(isapprox.(sol["solution"]["bus"]["1"]["vm"], [1.0189, 1.0313, 1.0313]; atol=2e-2))
            @test all(isapprox.(sol["solution"]["bus"]["1"]["va"]*180/pi, [30, -90, 150]; atol=3e-3))
            @test all(isapprox.(sol["solution"]["transformer"]["2"]["tap"], [1.01875, 1.03125, 1.03125]; atol=2e-2))
        end
        @testset "regcontrol_lpubfdiag" begin
            eng = parse_file("../test/data/opendss/IEEE13_RegControl.dss")
            math = transform_data_model(eng)
            sol = solve_mc_opf_oltc(math, LPUBFDiagPowerModel, ipopt_solver; make_si=false, solution_processors=[sol_data_model!])
            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"]), 0.405132; atol=7e-3)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"]), 0.266410; atol=1e-2)
            @test all(isapprox.(sol["solution"]["bus"]["1"]["vm"], [1.01545, 1.04077, 1.04216]; atol=3e-2))
            @test all(isapprox.(sol["solution"]["transformer"]["2"]["tap"], [1.01535, 1.04071, 1.04210]; atol=3e-2))
        end
        @testset "regcontrol_fbs" begin
            eng = parse_file("../test/data/opendss/IEEE13_RegControl.dss")
            math = transform_data_model(eng)
            sol = solve_mc_opf_oltc(math, FBSUBFPowerModel, ipopt_solver; make_si=false,solution_processors=[sol_data_model!])
            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"]), 0.405132; atol=7e-3)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"]), 0.266410; atol=1e-2)
            @test all(isapprox.(sol["solution"]["bus"]["1"]["vm"], [1.01545, 1.04077, 1.04216]; atol=3e-2))
            @test all(isapprox.(sol["solution"]["bus"]["1"]["va"]*180/pi, [29.9998, -90.0022, 149.9971]; atol=5e-3))
            @test all(isapprox.(sol["solution"]["transformer"]["2"]["tap"], [1.01535, 1.04071, 1.04210]; atol=3e-2))
        end
        @testset "regcontrol_fotp" begin
            eng = parse_file("../test/data/opendss/IEEE13_RegControl.dss")
            math = transform_data_model(eng)
            sol = solve_mc_opf_oltc(math, FOTPUPowerModel, ipopt_solver; make_si=false)
            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"]), 0.405132; atol=7e-3)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"]), 0.266410; atol=1e-2)
            @test all(isapprox.(sol["solution"]["bus"]["1"]["vm"], [1.01545, 1.04077, 1.04216]; atol=3e-2))
            @test all(isapprox.(sol["solution"]["bus"]["1"]["va"]*180/pi, [29.9998, -90.0022, 149.9971]; atol=5e-3))
            @test all(isapprox.(sol["solution"]["transformer"]["2"]["tap"], [1.01535, 1.04071, 1.04210]; atol=3e-2))
        end
        @testset "regcontrol_fotr" begin
            eng = parse_file("../test/data/opendss/IEEE13_RegControl.dss")
            math = transform_data_model(eng)
            sol = solve_mc_opf_oltc(math, FOTRUPowerModel, ipopt_solver; make_si=false,solution_processors=[sol_data_model!])
            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"]), 0.405132; atol=7e-3)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"]), 0.266410; atol=1e-2)
            @test all(isapprox.(sol["solution"]["bus"]["1"]["vm"], [1.01545, 1.04077, 1.04216]; atol=3e-2))
            @test all(isapprox.(sol["solution"]["bus"]["1"]["va"]*180/pi, [29.9998, -90.0022, 149.9971]; atol=5e-3))
            @test all(isapprox.(sol["solution"]["transformer"]["2"]["tap"], [1.01535, 1.04071, 1.04210]; atol=3e-2))
        end
        @testset "regcontrol_fotr" begin
            eng = parse_file("../test/data/opendss/IEEE13_RegControl.dss")
            math = transform_data_model(eng)
            sol = solve_mc_opf_oltc(math, FOTRUPowerModel, ipopt_solver; make_si=false,solution_processors=[sol_data_model!])
            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"]), 0.405132; atol=7e-3)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"]), 0.266410; atol=1e-2)
            @test all(isapprox.(sol["solution"]["bus"]["1"]["vm"], [1.01545, 1.04077, 1.04216]; atol=3e-2))
            @test all(isapprox.(sol["solution"]["bus"]["1"]["va"]*180/pi, [29.9998, -90.0022, 149.9971]; atol=5e-3))
            @test all(isapprox.(sol["solution"]["transformer"]["2"]["tap"], [1.01535, 1.04071, 1.04210]; atol=4e-2))
        end
    end
end
