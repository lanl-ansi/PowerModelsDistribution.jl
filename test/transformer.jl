# helper functions to access solutions by their OpenDSS names
bus_name2id(pmd_data, name) = [bus["index"] for (_,bus) in pmd_data["bus"] if haskey(bus, "name") && bus["name"]==name][1]
va(sol, pmd_data, name) = round.(PMD._wrap_to_pi(sol["solution"]["bus"][string(bus_name2id(pmd_data, name))]["va"][:])*180/pi; digits=1)
vm(sol, pmd_data, name) = sol["solution"]["bus"][string(bus_name2id(pmd_data, name))]["vm"]
# tests
@testset "transformer" begin
    @testset "two winding" begin
        @testset "ac_pf yy" begin
            file = "../test/data/opendss/ut_trans_2w_yy.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_tp_pf(pmd_data, ipopt_solver, multiconductor=true)
            @test norm(vm(sol, pmd_data, "3")-[0.87451, 0.8613, 0.85348], Inf) <= 1.5E-5
            @test norm(va(sol, pmd_data, "3")-[-0.1, -120.4, 119.8], Inf) <= 0.5E-2
        end
        @testset "ac_pf dy_lead" begin
            file = "../test/data/opendss/ut_trans_2w_dy_lead.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_tp_pf(pmd_data, ipopt_solver, multiconductor=true)
            @test norm(vm(sol, pmd_data, "3")-[0.87391, 0.86055, 0.85486], Inf) <= 1.5E-5
            @test norm(va(sol, pmd_data, "3")-[29.8, -90.4, 149.8], Inf) <= 0.5E-2
        end
        @testset "ac_pf dy_lag" begin
            file = "../test/data/opendss/ut_trans_2w_dy_lag.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_tp_pf(pmd_data, ipopt_solver, multiconductor=true)
            @test norm(vm(sol, pmd_data, "3")-[0.92092, 0.91012, 0.90059], Inf) <= 1.5E-5
            @test norm(va(sol, pmd_data, "3")-[-30.0, -150.4, 89.8], Inf) <= 0.5E-2
        end
        @testset "ac_opf_oltc" begin
            file = "../test/data/opendss/ut_trans_2w_yy_oltc.dss"
            pmd_data = PMD.parse_file(file)
            # free the taps
            pmd_data["trans"]["1"]["fixed"] = PMs.MultiConductorVector(zeros(Bool, 3))
            pmd_data["trans"]["2"]["fixed"] = PMs.MultiConductorVector(zeros(Bool, 3))
            pm = PMs.build_model(pmd_data, PMs.ACPPowerModel, PMD.post_tp_opf_oltc, ref_extensions=[PMD.ref_add_arcs_trans!], multiconductor=true)
            sol = PMs.optimize_model!(pm, ipopt_solver)
            # check that taps are set as to boost the voltage in the branches as much as possible;
            # this is trivially optimal if the voltage bounds are not binding
            # and without significant shunts (both branch and transformer)
            tap(i) = [JuMP.value(PMs.var(pm, pm.cnw, :cnd, c)[:tap][i]) for c in 1:3]
            @test norm(tap(1)-[0.95, 0.95, 0.95], Inf) <= 1E-4
            @test norm(tap(2)-[1.05, 1.05, 1.05], Inf) <= 1E-4
            # then check whether voltage is what OpenDSS expects for those values
            @test norm(vm(sol, pmd_data, "3")-[1.0352, 1.022, 1.0142], Inf) <= 1E-4
            @test norm(va(sol, pmd_data, "3")-[-0.1, -120.4, 119.8], Inf) <= 0.5E-2
        end
    end
    @testset "three winding" begin
        @testset "all  non-zero"  begin
            file = "../test/data/opendss/ut_trans_3w_dyy_1.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_tp_pf(pmd_data, ipopt_solver, multiconductor=true)
            @test norm(vm(sol, pmd_data, "3")-[0.9318, 0.88828, 0.88581], Inf) <= 1.5E-5
            @test norm(va(sol, pmd_data, "3")-[30.1, -90.7, 151.2], Inf) <= 0.5E-2
        end
        @testset "some non-zero" begin
            file = "../test/data/opendss/ut_trans_3w_dyy_2.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_tp_pf(pmd_data, ipopt_solver, multiconductor=true)
            #@test isapprox(vm(sol, pmd_data, "3"), [0.93876, 0.90227, 0.90454], atol=1E-4)
            @test norm(vm(sol, pmd_data, "3")-[0.93876, 0.90227, 0.90454], Inf) <= 1.5E-5
            @test norm(va(sol, pmd_data, "3")-[31.6, -88.8, 153.3], Inf) <= 0.5E-2
        end
        @testset "all zero" begin
            file = "../test/data/opendss/ut_trans_3w_dyy_3.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_tp_pf(pmd_data, ipopt_solver, multiconductor=true)
            @test norm(vm(sol, pmd_data, "3")-[0.97047, 0.93949, 0.946], Inf) <= 1.5E-5
            @test norm(va(sol, pmd_data, "3")-[30.6, -90.0, 151.9], Inf) <= 0.5E-2
        end
    end
    @testset "voltage base" begin
        # make sure that different voltage bases lead to the same solution
        # (after rescaling)
        file = "../test/data/opendss/ut_trans_3w_dyy_basetest.dss"
        pmd1 = PMD.parse_file(file)
        pmd2 = deepcopy(pmd1)
        PMD._adjust_base!(pmd2, start_at_first_tr_prim=false)
        scale_2to1 = pmd2["bus"]["3"]["base_kv"]/pmd1["bus"]["3"]["base_kv"]
        sol1 = PMD.run_ac_tp_pf(pmd1, ipopt_solver, multiconductor=true)
        sol2 = PMD.run_ac_tp_pf(pmd2, ipopt_solver, multiconductor=true)
        for bus_id_str in keys(sol1["solution"]["bus"])
            vm1 = sol1["solution"]["bus"][bus_id_str]["vm"]
            vm2 = sol2["solution"]["bus"][bus_id_str]["vm"]
            va1 = PMD._wrap_to_pi(sol1["solution"]["bus"][bus_id_str]["va"][:])
            va2 = PMD._wrap_to_pi(sol2["solution"]["bus"][bus_id_str]["va"][:])
            @test norm(vm1-vm2*scale_2to1, Inf) <= 1E-6
            @test norm(va1-va2, Inf) <= 1E-3
        end
    end
end
