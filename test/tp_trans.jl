# helper functions
bus_name2id(tppm_data, name) = [bus["index"] for (_,bus) in tppm_data["bus"] if haskey(bus, "name") && bus["name"]==name][1]
va(sol, tppm_data, name) = round.(TPPMs.wraptopi(sol["solution"]["bus"][string(bus_name2id(tppm_data, name))]["va"])*180/pi; digits=1)
vm(sol, tppm_data, name) = sol["solution"]["bus"][string(bus_name2id(tppm_data, name))]["vm"]
# tests
@testset "transformer" begin
    @testset "two winding" begin
        @testset "ac_pf yy" begin
            file = "../test/data/opendss/ut_trans_2w_yy.dss"
            tppm_data = TPPMs.parse_file(file)
            sol = TPPMs.run_ac_tp_pf(tppm_data, ipopt_solver, multiconductor=true)
            @test norm(vm(sol, tppm_data, "3")-[0.87451, 0.8613, 0.85348], Inf) <= 1.5E-5
            @test norm(va(sol, tppm_data, "3")-[-0.1, -120.4, 119.8], Inf) <= 0.5E-2
        end
        @testset "ac_pf dy_lead" begin
            file = "../test/data/opendss/ut_trans_2w_dy_lead.dss"
            tppm_data = TPPMs.parse_file(file)
            sol = TPPMs.run_ac_tp_pf(tppm_data, ipopt_solver, multiconductor=true)
            @test norm(vm(sol, tppm_data, "3")-[0.87391, 0.86055, 0.85486], Inf) <= 1.5E-5
            @test norm(va(sol, tppm_data, "3")-[29.8, -90.4, 149.8], Inf) <= 0.5E-2
        end
        @testset "ac_pf dy_lag" begin
            file = "../test/data/opendss/ut_trans_2w_dy_lag.dss"
            tppm_data = TPPMs.parse_file(file)
            sol = TPPMs.run_ac_tp_pf(tppm_data, ipopt_solver, multiconductor=true)
            @test norm(vm(sol, tppm_data, "3")-[0.92092, 0.91012, 0.90059], Inf) <= 1.5E-5
            @test norm(va(sol, tppm_data, "3")-[-30.0, -150.4, 89.8], Inf) <= 0.5E-2
        end
        @testset "ac_opf oltc" begin
            file = "../test/data/opendss/ut_trans_2w_yy_oltc.dss"
            tppm_data = TPPMs.parse_file(file)
            # free the taps
            tppm_data["trans"]["1"]["tapfix"] = MultiConductorVector(zeros(Bool, 3))
            tppm_data["trans"]["2"]["tapfix"] = MultiConductorVector(zeros(Bool, 3))
            pm = PMs.build_generic_model(tppm_data, PMs.ACPPowerModel, TPPMs.post_tp_opf, multiconductor=true)
            sol = PMs.solve_generic_model(pm, ipopt_solver)
            # check that taps are set as to boost the voltage in the branches as much as possible;
            # this is trivially optimal if the voltage bounds are not binding
            # and without significant shunts (both branch and transformer)
            tap(i) = [JuMP.getvalue(var(pm, pm.cnw, :cnd, c)[:tap][i]) for c in 1:3]
            @test norm(tap(1)-[0.95, 0.95, 0.95], Inf) <= 1E-4
            @test norm(tap(2)-[1.05, 1.05, 1.05], Inf) <= 1E-4
            # then check whether voltage is what OpenDSS expects for those values
            @test norm(vm(sol, tppm_data, "3")-[1.0352, 1.022, 1.0142], Inf) <= 1E-4
            @test norm(va(sol, tppm_data, "3")-[-0.1, -120.4, 119.8], Inf) <= 0.5E-2
        end
    end
    @testset "three winding" begin
        @testset "all  non-zero"  begin
            file = "../test/data/opendss/ut_trans_3w_dyy_1.dss"
            tppm_data = TPPMs.parse_file(file)
            sol = TPPMs.run_ac_tp_pf(tppm_data, ipopt_solver, multiconductor=true)
            @test norm(vm(sol, tppm_data, "3")-[0.9318, 0.88828, 0.88581], Inf) <= 1.5E-5
            @test norm(va(sol, tppm_data, "3")-[30.1, -90.7, 151.2], Inf) <= 0.5E-2
        end
        @testset "some non-zero" begin
            file = "../test/data/opendss/ut_trans_3w_dyy_2.dss"
            tppm_data = TPPMs.parse_file(file)
            sol = TPPMs.run_ac_tp_pf(tppm_data, ipopt_solver, multiconductor=true)
            #@test isapprox(vm(sol, tppm_data, "3"), [0.93876, 0.90227, 0.90454], atol=1E-4)
            @test norm(vm(sol, tppm_data, "3")-[0.93876, 0.90227, 0.90454], Inf) <= 1.5E-5
            @test norm(va(sol, tppm_data, "3")-[31.6, -88.8, 153.3], Inf) <= 0.5E-2
        end
        @testset "all zero" begin
            file = "../test/data/opendss/ut_trans_3w_dyy_3.dss"
            tppm_data = TPPMs.parse_file(file)
            sol = TPPMs.run_ac_tp_pf(tppm_data, ipopt_solver, multiconductor=true)
            @test norm(vm(sol, tppm_data, "3")-[0.97047, 0.93949, 0.946], Inf) <= 1.5E-5
            @test norm(va(sol, tppm_data, "3")-[30.6, -90.0, 151.9], Inf) <= 0.5E-2
        end
    end
end
