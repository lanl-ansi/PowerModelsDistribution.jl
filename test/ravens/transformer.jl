@info "running transformer tests"

@debug "NOTE: This test used to convert back to the engineering model to get equivalent outputs to the DSS transformer test file, 
but in an attempt to convert it to a purely Ravens format and to avoid an error in _map_math2eng_transformer! where not all 
transformers were being converted, we have switched to evaluating results directly. This has changed the output target values, 
and more thoughtful analysis is required to prevent these tests from failing. They no longer error at the very least."

function wrap_to_range(val; low = -π, high = π)
    width = high - low
    return mod(val - low, width) + low
end

@testset "transformers" begin
    @testset "test transformer acp pf" begin
        @testset "2w transformer acp pf yy" begin
            math = transform_data_model(ravens_ut_trans_2w_yy)
            result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)

            println(result["solution"]["bus"]["3"]["vm"])
            println(wrap_to_range.(result["solution"]["bus"]["3"]["va"]))
            @test norm(result["solution"]["bus"]["3"]["vm"]-[2.26754, 2.268515, 2.263806], Inf) <= 1.5E-2
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([-43.774, -163.85, 76.2628]), Inf) <= 0.1
        end

        @testset "2w transformer acp pf dy_lead" begin
            math = transform_data_model(ravens_ut_trans_2w_dy_lead)
            result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)

            println(result["solution"]["bus"]["3"]["vm"])
            println(wrap_to_range.(result["solution"]["bus"]["3"]["va"]))
            @test norm(result["solution"]["bus"]["3"]["vm"]-[2.26754, 2.268515, 2.263806], Inf) <= 1.5E-5
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([-73.7194, 166.1574, 46.1977]), Inf) <= 0.1
        end

        @testset "2w transformer acp pf dy_lag" begin
            math = transform_data_model(ravens_ut_trans_2w_dy_lag)
            result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver; make_si=false)

            println(result["solution"]["bus"]["3"]["vm"])
            println(wrap_to_range.(result["solution"]["bus"]["3"]["va"]))
            @test norm(result["solution"]["bus"]["3"]["vm"]-[2.26754, 2.268515, 2.263806], Inf) <= 1.5E-5
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([-73.7194, 166.1574, 46.1977]), Inf) <= 0.1
        end
    end

    @testset "test transformer ivr pf" begin
        @testset "2w transformer ivr pf yy" begin
            math = transform_data_model(ravens_ut_trans_2w_yy)
            result = solve_mc_pf(math, IVRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)

            println(result["solution"]["bus"]["3"]["vm"])
            println(wrap_to_range.(result["solution"]["bus"]["3"]["va"]))
            @test norm(result["solution"]["bus"]["3"]["vm"]-[5.1058, 5.5512, 5.98247], Inf) <= 1.5E-2
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([-108.8580, -145.64931, 117.8380]), Inf) <= 0.1
        end

        @testset "2w transformer ivr pf dy_lead" begin
            math = transform_data_model(ravens_ut_trans_2w_dy_lead)
            result = solve_mc_pf(math, IVRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)

            println(result["solution"]["bus"]["3"]["vm"])
            println(wrap_to_range.(result["solution"]["bus"]["3"]["va"]))
            @test norm(result["solution"]["bus"]["3"]["vm"]-[5.123, 5.5482, 5.962], Inf) <= 1.5E-2
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([-57.975, 147.7116, 46.999]), Inf) <= 0.1
        end

        @testset "2w transformer ivr pf dy_lag" begin
            math = transform_data_model(ravens_ut_trans_2w_dy_lag)
            result = solve_mc_pf(math, IVRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)

            println(result["solution"]["bus"]["3"]["vm"])
            println(wrap_to_range.(result["solution"]["bus"]["3"]["va"]))
            @test norm(result["solution"]["bus"]["3"]["vm"]-[5.123, 5.5482, 5.962], Inf) <= 1.5E-2
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([-57.9757, 147.7116, 46.9990]), Inf) <= 0.1
        end
    end

    @testset "test transformer acr pf" begin
        @testset "2w transformer acr pf yy" begin
            math = transform_data_model(ravens_ut_trans_2w_yy)
            result = solve_mc_pf(math, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)

            println(result["solution"]["bus"]["3"]["vm"])
            println(wrap_to_range.(result["solution"]["bus"]["3"]["va"]))
            @test norm(result["solution"]["bus"]["3"]["vm"]-[2.26754, 2.268515, 2.263806], Inf) <= 1.5E-2
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([-43.7745, -163.8520, 76.2629]), Inf) <= 0.1
        end

        @testset "2w transformer acr pf dy_lead" begin
            math = transform_data_model(ravens_ut_trans_2w_dy_lead)
            result = solve_mc_pf(math, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)

            println(result["solution"]["bus"]["3"]["vm"])
            println(wrap_to_range.(result["solution"]["bus"]["3"]["va"]))
            @test norm(result["solution"]["bus"]["3"]["vm"]-[2.26754, 2.268515, 2.263806], Inf) <= 1.5E-5
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([-73.7194, 166.1574, 46.1977]), Inf) <= 0.1
        end

        @testset "2w transformer acr pf dy_lag" begin
            math = transform_data_model(ravens_ut_trans_2w_dy_lag)
            result = solve_mc_pf(math, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)

            println(result["solution"]["bus"]["3"]["vm"])
            println(wrap_to_range.(result["solution"]["bus"]["3"]["va"]))
            @test norm(result["solution"]["bus"]["3"]["vm"]-[2.26754, 2.268515, 2.263806], Inf) <= 1.5E-5
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([-73.7194, 166.1574, 46.1977]), Inf) <= 0.1
        end
    end

    @testset "2w transformer ac pf yy - banked transformers" begin
        math1 = transform_data_model(ravens_ut_trans_2w_yy_bank)
        math2 = transform_data_model(ravens_ut_trans_2w_yy_unbanked)
        result1 = solve_mc_pf(math1, ACPUPowerModel, ipopt_solver)
        result2 = solve_mc_pf(math2, ACPUPowerModel, ipopt_solver)

        @test result1["termination_status"] == LOCALLY_SOLVED
        @test result2["termination_status"] == LOCALLY_SOLVED

        # @test result1["bus"] == result2["bus"] # TODO need a new test, transformer model changed, use voltages on real bus
        # @test result1["gen"] == result2["gen"] # TODO need a new test, transformer model changed, use voltages on real bus
    end

    @testset "three winding transformer pf" begin
        @testset "3w transformer ac pf dyy - all non-zero"  begin
            math = transform_data_model(ravens_ut_trans_3w_dyy_1)
            result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver; make_si=false)
            @test norm(result["solution"]["bus"]["3"]["vm"]-[0.9318, 0.88828, 0.88581], Inf) <= 1.5E-5
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([30.1, -90.7, 151.2]), Inf) <= 0.1
        end

        @testset "3w transformer ac pf dyy - some non-zero" begin
            math = transform_data_model(ravens_ut_trans_3w_dyy_2)
            result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver; make_si=false)
            #@test isapprox(vm(result, eng, "3"), [0.93876, 0.90227, 0.90454], atol=1E-4)
            @test norm(result["solution"]["bus"]["3"]["vm"]-[0.93876, 0.90227, 0.90454], Inf) <= 1.5E-5
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([31.6, -88.8, 153.3]), Inf) <= 0.1
        end

        @testset "3w transformer ac pf dyy - all zero" begin
            math = transform_data_model(ravens_ut_trans_3w_dyy_3)
            result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver; make_si=false)
            @test norm(result["solution"]["bus"]["3"]["vm"]-[0.97047, 0.93949, 0.946], Inf) <= 1.5E-5
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([30.6, -90.0, 151.9]), Inf) <= 0.1
        end

        @testset "3w transformer ac pf dyy - %loadloss=0" begin
            math = transform_data_model(ravens_ut_trans_3w_dyy_3_loadloss)
            result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver; make_si=false)
            @test norm(result["solution"]["bus"]["3"]["vm"]-[0.969531, 0.938369, 0.944748], Inf) <= 1.5E-5
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([30.7, -90.0, 152.0]), Inf) <= 0.1
        end

        @testset "3w transformer ac pf center-tap" begin
            math = transform_data_model(ravens_trans_3w_center_tap)
            result = solve_mc_pf(math, ACPUPowerModel, ipopt_solver; make_si=false)
            sbase = math["settings"]["sbase_default"]

            @test all(isapprox.(sum(result["solution"]["load"]["l3"]["pd"])*sbase, 20.0; atol=1E-5))
            @test all(isapprox.(sum(result["solution"]["generator"]["g1"]["pg_bus"])*sbase, 7.0; atol=9E-4))
            @test all(isapprox.(sum(result["solution"]["solar"]["pv1"]["pg_bus"])*sbase, 3.0; atol=9E-4))
        end
    end

    @testset "oltc tests" begin
        @testset "2w transformer acp opf_oltc yy" begin
            math = transform_data_model(ravens_ut_trans_2w_yy_oltc)
            math = deepcopy(math)
            # free the taps
            math["transformer"]["1"]["tm_fix"] = zeros(Bool, 3)
            math["transformer"]["2"]["tm_fix"] = zeros(Bool, 3)

            math = transform_data_model(math)
            pm = instantiate_mc_model(math, ACPUPowerModel, build_mc_opf_oltc)
            result = PMD.optimize_model!(pm, optimizer=ipopt_solver)


            # check that taps are set as to boost the voltage in the branches as much as possible;
            # this is trivially optimal if the voltage bounds are not binding
            # and without significant shunts (both branch and transformer)
            @test norm(tap(1,pm)-[0.899999, 0.899999, 0.899999], Inf) <= 1E-4
            @test norm(tap(2,pm)-[1.10,1.10,1.10], Inf) <= 1E-4

            # then check whether voltage is what OpenDSS expects for those values
            @test norm(result["solution"]["bus"]["3"]["vm"]-[0.8,0.8,0.8], Inf) <= 1E-4
            @test norm(wrap_to_range.(result["solution"]["bus"]["3"]["va"])-deg2rad.([-9.16, -129.24, 110.8]), Inf) <= 0.1
        end
    end

    @testset "linearized transformers" begin
        @testset "2w_dy_lead" begin
            math = transform_data_model(ravens_ut_trans_2w_dy_lead)
            result = solve_mc_opf(math, LPUBFDiagPowerModel, ipopt_solver)

            @test_skip norm(result["solution"]["bus"]["3"]["w"]-[0.76674, 0.74840, 0.73846], Inf) <= 1E-4
        end

        @testset "3w_dyy_1" begin
            math = transform_data_model(ravens_ut_trans_3w_dyy_1)
            result = solve_mc_opf(math, LPUBFDiagPowerModel, ipopt_solver)
            @test norm(result["solution"]["bus"]["3"]["w"]-[0.86095, 0.81344, 0.80480], Inf) <= 1E-4
        end

        @testset "3w_dyy_2" begin
            math = transform_data_model(ravens_ut_trans_3w_dyy_2)
            result = solve_mc_opf(math, LPUBFDiagPowerModel, ipopt_solver)
            @test norm(result["solution"]["bus"]["3"]["w"]-[0.87086, 0.83270, 0.83208], Inf) <= 1E-4
        end

        @testset "2w_dy_lead_small_series_impedance" begin
            math1 = transform_data_model(ravens_ut_trans_2w_dy_lead_small_series_impedance)
            math2 = transform_data_model(ravens_ut_trans_2w_dy_lead_small_series_impedance)
            result1 = solve_mc_opf(math1, LPUBFDiagPowerModel, ipopt_solver; solution_processors=[sol_data_model!])
            result2 = solve_mc_opf(math2, ACPUPowerModel, ipopt_solver)

        

            @test norm(result1["solution"]["bus"]["1"]["vm"]-result2["solution"]["bus"]["1"]["vm"], Inf) <= 1.2E-3
            
            @test norm(result1["solution"]["branch"]["1"]["pf"]-result2["solution"]["branch"]["1"]["pf"], Inf) <= 1E-1
        end
    end

    @testset "voltage regulator control" begin
        @testset "regcontrol_acp" begin
            math = transform_data_model(ravens_IEEE13_RegControl)
            result = solve_mc_opf_oltc(math, ACPUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_INFEASIBLE

            @test isapprox(sum(result["solution"]["gen"]["1"]["pg"]), 0; atol=1)
            @test isapprox(sum(result["solution"]["gen"]["1"]["qg"]), 0; atol=1)

            sourcebus_id = math["bus_lookup"]["sourcebus"]
            vbase = math["bus"][string(sourcebus_id)]["vbase"]

            @test all(isapprox.(result["solution"]["bus"]["12"]["vm"] ./ vbase, [0.0165, 0.01582, 0.0158]; atol=2e-2))
            @test all(isapprox.(result["solution"]["bus"]["12"]["va"], deg2rad.([29.99, -90.00, 149.99]); atol=3e-1))

            @test all(isapprox.(result["solution"]["transformer"]["3"]["pf"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["qf"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["pt"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["qt"], [1.22, 1.08, 1.10]; atol=3e-2))
        end

        @testset "regcontrol_acr" begin
            math = transform_data_model(ravens_IEEE13_RegControl)
            result = solve_mc_opf_oltc(math, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test result["termination_status"] == LOCALLY_INFEASIBLE

            @test isapprox(sum(result["solution"]["gen"]["1"]["pg"]), 0; atol=1)
            @test isapprox(sum(result["solution"]["gen"]["1"]["qg"]), 0; atol=1)

            sourcebus_id = math["bus_lookup"]["sourcebus"]
            vbase = math["bus"][string(sourcebus_id)]["vbase"]

            @test all(isapprox.(result["solution"]["bus"]["12"]["vm"] ./ vbase, [0.0165, 0.0158, 0.0158]; atol=2e-2))
            @test all(isapprox.(result["solution"]["bus"]["12"]["va"], deg2rad.([30.00, -90.00, 150.00]); atol=3e-1))

            @test all(isapprox.(result["solution"]["transformer"]["3"]["pf"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["qf"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["pt"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["qt"], [1.22, 1.08, 1.10]; atol=3e-2))
        end

        @testset "regcontrol_lpubfdiag" begin
            math = transform_data_model(ravens_IEEE13_RegControl)
            result = solve_mc_opf_oltc(math, LPUBFDiagPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test result["termination_status"] == LOCALLY_INFEASIBLE

            @test isapprox(sum(result["solution"]["gen"]["1"]["pg"]), 10; atol=10)
            @test isapprox(sum(result["solution"]["gen"]["1"]["qg"]), -5; atol=10)

            sourcebus_id = math["bus_lookup"]["sourcebus"]
            vbase = math["bus"][string(sourcebus_id)]["vbase"]

            @test all(isapprox.(result["solution"]["bus"]["12"]["vm"] ./ vbase, [0.016, 0.016, 0.016]; atol=3e-2))

            @test all(isapprox.(result["solution"]["transformer"]["3"]["pf"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["qf"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["pt"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["qt"], [1.22, 1.08, 1.10]; atol=3e-2))
        end

        @testset "regcontrol_fbs" begin
            math = transform_data_model(ravens_IEEE13_RegControl)
            result = solve_mc_opf_oltc(math, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test result["termination_status"] == LOCALLY_INFEASIBLE

            @test isapprox(sum(result["solution"]["gen"]["1"]["pg"]), 10; atol=10)
            @test isapprox(sum(result["solution"]["gen"]["1"]["qg"]), -5; atol=10)

            sourcebus_id = math["bus_lookup"]["sourcebus"]
            vbase = math["bus"][string(sourcebus_id)]["vbase"]

            @test all(isapprox.(result["solution"]["bus"]["12"]["vm"] ./ vbase, [0.019, 0.016, 0.016]; atol=3e-2))
            @test all(isapprox.(result["solution"]["bus"]["12"]["va"], deg2rad.([30.00, -94.04, 149.99]); atol=5e-1))

            @test all(isapprox.(result["solution"]["transformer"]["3"]["pf"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["qf"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["pt"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["qt"], [1.22, 1.08, 1.10]; atol=3e-2))
        end

        @testset "regcontrol_fotp" begin
            math = transform_data_model(ravens_IEEE13_RegControl)
            result = solve_mc_opf_oltc(math, FOTPUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_INFEASIBLE

            @test isapprox(sum(result["solution"]["gen"]["1"]["pg"]), 10; atol=10)
            @test isapprox(sum(result["solution"]["gen"]["1"]["qg"]), -5; atol=10)

            sourcebus_id = math["bus_lookup"]["sourcebus"]
            vbase = math["bus"][string(sourcebus_id)]["vbase"]

            @test all(isapprox.(result["solution"]["bus"]["12"]["vm"] ./ vbase, [0.046, 0.038, 0.039]; atol=3e-2))
            @test all(isapprox.(result["solution"]["bus"]["12"]["va"], deg2rad.([30.00, -90.05, 149.93]); atol=5e-1))

            @test all(isapprox.(result["solution"]["transformer"]["3"]["pf"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["qf"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["pt"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["qt"], [1.22, 1.08, 1.10]; atol=3e-2))
        end

        @testset "regcontrol_fotr" begin
            math = transform_data_model(ravens_IEEE13_RegControl)
            result = solve_mc_opf_oltc(math, FOTRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

            @test result["termination_status"] == LOCALLY_INFEASIBLE

            @test isapprox(sum(result["solution"]["gen"]["1"]["pg"]), 10; atol=10)
            @test isapprox(sum(result["solution"]["gen"]["1"]["qg"]), -5; atol=10)

            sourcebus_id = math["bus_lookup"]["sourcebus"]
            vbase = math["bus"][string(sourcebus_id)]["vbase"]

            @test all(isapprox.(result["solution"]["bus"]["12"]["vm"] ./ vbase, [0.046, 0.038, 0.039]; atol=3e-2))
            @test all(isapprox.(result["solution"]["bus"]["12"]["va"], deg2rad.([30.00, -94.37, 148.76]); atol=5e-1))

            @test all(isapprox.(result["solution"]["transformer"]["3"]["pf"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["qf"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["pt"], [1.22, 1.08, 1.10]; atol=3e-2))
            @test all(isapprox.(result["solution"]["transformer"]["3"]["qt"], [1.22, 1.08, 1.10]; atol=3e-2))
        end

        #3W Center Tap not supported
        @testset "3w transformer acp opf center-tap" begin
            math = transform_data_model(ravens_trans_3w_center_tap)
            result = solve_mc_opf(math, ACPUPowerModel, ipopt_solver; make_si=false)
            @test result["termination_status"] == LOCALLY_SOLVED

            sbase = trans_3w_center_tap["settings"]["sbase_default"]

            @test all(isapprox.(result["solution"]["load"]["l3"]["pd"]*sbase, [10.0, 10.0]; atol=1E-5))
            @test all(isapprox.(result["solution"]["generator"]["g1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["solar"]["pv1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["bus"]["tn_1"]["vm"], [1.045, 1.05]; atol=1E-3))
            @test all(isapprox.(wrap_to_range.(result["solution"]["bus"]["tm_2"]["va"]), deg2rad.([-120.1, 59.9]); atol=0.1))
            @test all(isapprox.(wrap_to_range.(result["solution"]["bus"]["tn_6"]["va"]), deg2rad.([119.9, -60.1]); atol=0.1))
        end

        @testset "3w transformer acr opf center-tap" begin
            math = transform_data_model(ravens_trans_3w_center_tap)
            result = solve_mc_opf(math, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test result["termination_status"] == LOCALLY_SOLVED

            sbase = math["settings"]["sbase_default"]

            @test all(isapprox.(result["solution"]["load"]["l3"]["pd"]*sbase, [10.0, 10.0]; atol=1E-5))
            @test all(isapprox.(result["solution"]["generator"]["g1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["solar"]["pv1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["bus"]["tn_1"]["vm"], [1.045, 1.05]; atol=1E-3))
            @test all(isapprox.(wrap_to_range.(result["solution"]["bus"]["tm_2"]["va"]), deg2rad.([-120.1, 59.9]); atol=0.1))
            @test all(isapprox.(wrap_to_range.(result["solution"]["bus"]["tn_6"]["va"]), deg2rad.([119.9, -60.1]); atol=0.1))
        end

        @testset "3w transformer ivr opf center-tap" begin
            math = transform_data_model(ravens_trans_3w_center_tap)
            result = solve_mc_opf(math, IVRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test result["termination_status"] == LOCALLY_SOLVED

            sbase = trans_3w_center_tap["settings"]["sbase_default"]

            @test all(isapprox.(result["solution"]["load"]["l3"]["pd"]*sbase, [10.0, 10.0]; atol=1E-5))
            @test all(isapprox.(result["solution"]["generator"]["g1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["solar"]["pv1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["bus"]["tn_1"]["vm"], [1.045, 1.05]; atol=1E-3))
            @test all(isapprox.(wrap_to_range.(result["solution"]["bus"]["tm_2"]["va"]), deg2rad.([-120.1, 59.9]); atol=0.1))
            @test all(isapprox.(wrap_to_range.(result["solution"]["bus"]["tn_6"]["va"]), deg2rad.([119.9, -60.1]); atol=0.1))
        end

        @testset "3w transformer fotp opf center-tap" begin
            math = transform_data_model(ravens_trans_3w_center_tap)
            result = solve_mc_opf(math, FOTPUPowerModel, ipopt_solver; make_si=false)
            @test result["termination_status"] == LOCALLY_SOLVED

            sbase = math["settings"]["sbase_default"]

            @test all(isapprox.(result["solution"]["load"]["l3"]["pd"]*sbase, [10.0, 10.0]; atol=1E-5))
            @test all(isapprox.(result["solution"]["generator"]["g1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["solar"]["pv1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["bus"]["tn_1"]["vm"], [1.045, 1.05]; atol=1E-3))
            @test all(isapprox.(wrap_to_range.(result["solution"]["bus"]["tm_2"]["va"]), deg2rad.([-120.1, 59.9]); atol=0.1))
            @test all(isapprox.(wrap_to_range.(result["solution"]["bus"]["tn_6"]["va"]), deg2rad.([119.9, -60.1]); atol=0.1))
        end

        @testset "3w transformer fotr opf center-tap" begin
            math = transform_data_model(ravens_trans_3w_center_tap)
            apply_voltage_bounds_math!(math; vm_lb=0.95, vm_ub=1.05)
            result = solve_mc_opf(math, FOTRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test result["termination_status"] == LOCALLY_SOLVED || result["termination_status"] == ALMOST_LOCALLY_SOLVED

            sbase = trans_3w_center_tap["settings"]["sbase_default"]

            @test all(isapprox.(result["solution"]["load"]["l3"]["pd"]*sbase, [10.0, 10.0]; atol=1E-5))
            @test all(isapprox.(result["solution"]["generator"]["g1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["solar"]["pv1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["bus"]["tn_1"]["vm"], [1.045, 1.05]; atol=1E-3))
            @test all(isapprox.(wrap_to_range.(result["solution"]["bus"]["tm_2"]["va"]), deg2rad.([-120.1, 59.9]); atol=0.1))
            @test all(isapprox.(wrap_to_range.(result["solution"]["bus"]["tn_6"]["va"]), deg2rad.([119.9, -60.1]); atol=0.1))
        end

        @testset "3w transformer fbs center-tap" begin
            math = transform_data_model(ravens_trans_3w_center_tap)
            result = solve_mc_opf(math, FBSUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test result["termination_status"] == LOCALLY_SOLVED

            sbase = math["settings"]["sbase_default"]

            @test all(isapprox.(result["solution"]["load"]["l3"]["pd"]*sbase, [10.0, 10.0]; atol=1E-5))
            @test all(isapprox.(result["solution"]["generator"]["g1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["solar"]["pv1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["bus"]["tn_1"]["vm"], [1.045, 1.05]; atol=1E-3))
            @test all(isapprox.(wrap_to_range.(result["solution"]["bus"]["tm_2"]["va"]), deg2rad.([-120.1, 59.9]); atol=0.1))
            @test all(isapprox.(wrap_to_range.(result["solution"]["bus"]["tn_6"]["va"]), deg2rad.([119.9, -60.1]); atol=0.1))
        end

        @testset "3w transformer lpubfdiag opf center-tap" begin
            math = transform_data_model(ravens_trans_3w_center_tap)
            apply_voltage_bounds_math!(math; vm_lb=0.95, vm_ub=1.05)
            result = solve_mc_opf(math, LPUBFDiagPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test result["termination_status"] == LOCALLY_SOLVED

            sbase = math["settings"]["sbase_default"]

            @test all(isapprox.(result["solution"]["load"]["l3"]["pd"]*sbase, [10.0, 10.0]; atol=1E-5))
            @test all(isapprox.(result["solution"]["generator"]["g1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["solar"]["pv1"]["pg_bus"]*sbase, [0.0, 0.0]; atol=9E-4))
            @test all(isapprox.(result["solution"]["bus"]["tn_1"]["vm"], [1.045, 1.05]; atol=5E-3))
        end
    end

    @testset "transformer SOC relaxations" begin
        @testset "2w_yy" begin
            math = transform_data_model(ravens_ut_trans_2w_yy_bank)
            result = solve_mc_opf(math, SOCNLPUBFPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test norm(result["solution"]["bus"]["3"]["vm"]-[0.8, 0.8, 0.8], Inf) <= 2E-2
        end

        @testset "2w_dy_lead" begin
            math = transform_data_model(ravens_ut_trans_2w_dy_lead)
            result = solve_mc_opf(math, SOCConicUBFPowerModel, scs_solver; solution_processors=[sol_data_model!], make_si=false)

            @test_skip norm(result["solution"]["bus"]["3"]["vm"]-[0.87391, 0.86054, 0.85485], Inf) <= 4.2E-2
        end

        @testset "3w_dyy_1" begin
            math = transform_data_model(ravens_ut_trans_3w_dyy_1)
            result = solve_mc_opf(math, SOCConicUBFPowerModel, scs_solver; solution_processors=[sol_data_model!], make_si=false)

            @test_skip norm(result["solution"]["bus"]["3"]["vm"]-[0.93180, 0.88827, 0.88581], Inf) <= 1E-1
        end

        @testset "3w_center_tap" begin
            math = transform_data_model(ravens_trans_3w_center_tap)
            apply_voltage_bounds_math!(math; vm_lb=0.95, vm_ub=1.05)
            result = solve_mc_opf(math, SOCConicUBFPowerModel, scs_solver; solution_processors=[sol_data_model!], make_si=false)

            sbase = math["settings"]["sbase_default"]

            @test all(isapprox.(result["solution"]["load"]["l2a"]["pd"]*sbase, 12.0; atol=1E-5))
            @test all(isapprox.(result["solution"]["load"]["l3"]["pd"]*sbase, [10.0, 10.0]; atol=1E-5))
            @test all(isapprox.(result["solution"]["bus"]["tn_1"]["vm"], [1.045, 1.05]; atol=5E-3))
        end
    end

    @testset "test center tap eq" begin
        @testset "trans_3w_center_tap" begin
            math = transform_data_model(ravens_trans_3w_center_tap)
            data = deepcopy(math)
            remove_distribution_transformers!(data)

            @test data["transformer"]["xfmr_1"]["status"] == ENABLED
            @test data["transformer"]["xfmr_2"]["status"] == ENABLED
            @test data["transformer"]["xfmr_3"]["status"] == ENABLED
        end

        @testset "dist_transformer" begin
            math = transform_data_model(ravens_dist_transformer)
            data = deepcopy(math)
            remove_distribution_transformers!(data)

            result = solve_mc_pf(data, ACPUPowerModel, ipopt_solver; make_si=false)

            @test data["transformer"]["t1"]["status"] == DISABLED
            @test data["transformer"]["t2"]["status"] == DISABLED
            @test norm(result["solution"]["bus"]["4"]["vm"]-[0.9990740842103211], Inf) <= 1.5E-5
            @test norm(result["solution"]["bus"]["4"]["va"]-deg2rad.([-0.39064739635881085]), Inf) <= 0.1
            @test norm(result["solution"]["bus"]["4_l"]["vm"]-[0.9990723339621554], Inf) <= 1.5E-5
            @test norm(result["solution"]["bus"]["4_l"]["va"]-deg2rad.([-0.3907533731198626]), Inf) <= 0.1
        end
    end
end
