@info "running transformer tests"

@testset "transformers" begin

    @testset "test transformer acp pf" begin
        @testset "2w transformer acp pf yy" begin
            file = "../test/data/opendss/ut_trans_2w_yy.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_mc_pf(pmd_data, ipopt_solver, multiconductor=true, solution_processors=[PMD.sol_polar_voltage!])
            solution_identify!(sol["solution"], pmd_data)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.87451, 0.8613, 0.85348], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-deg2rad.([-0.1, -120.4, 119.8]), Inf) <= 0.1
        end

        @testset "2w transformer acp pf dy_lead" begin
            file = "../test/data/opendss/ut_trans_2w_dy_lead.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_mc_pf(pmd_data, ipopt_solver, multiconductor=true, solution_processors=[PMD.sol_polar_voltage!])
            solution_identify!(sol["solution"], pmd_data)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.87391, 0.86055, 0.85486], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-deg2rad.([29.8, -90.4, 149.8]), Inf) <= 0.1
        end

        @testset "2w transformer acp pf dy_lag" begin
            file = "../test/data/opendss/ut_trans_2w_dy_lag.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_mc_pf(pmd_data, ipopt_solver, multiconductor=true)
            solution_identify!(sol["solution"], pmd_data)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.92092, 0.91012, 0.90059], Inf) <= 1.5E-5
            @test norm(PMD._wrap_to_pi(sol["solution"]["bus"]["3"]["va"])-deg2rad.([-30.0, -150.4, 89.8]), Inf) <= 0.1
        end
    end

    @testset "test transformer ivr pf" begin
        @testset "2w transformer ivr pf yy" begin
            file = "../test/data/opendss/ut_trans_2w_yy.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_mc_pf_iv(pmd_data, PMs.IVRPowerModel, ipopt_solver, solution_processors=[PMD.sol_polar_voltage!])
            solution_identify!(sol["solution"], pmd_data)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.87451, 0.8613, 0.85348], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-deg2rad.([-0.1, -120.4, 119.8]), Inf) <= 0.1
        end

        @testset "2w transformer ivr pf dy_lead" begin
            file = "../test/data/opendss/ut_trans_2w_dy_lead.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_mc_pf_iv(pmd_data, PMs.IVRPowerModel, ipopt_solver, solution_processors=[PMD.sol_polar_voltage!])
            solution_identify!(sol["solution"], pmd_data)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.87391, 0.86055, 0.85486], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-deg2rad.([29.8, -90.4, 149.8]), Inf) <= 0.1
        end

        @testset "2w transformer ivr pf dy_lag" begin
            file = "../test/data/opendss/ut_trans_2w_dy_lag.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_mc_pf_iv(pmd_data, PMs.IVRPowerModel, ipopt_solver, solution_processors=[PMD.sol_polar_voltage!])
            solution_identify!(sol["solution"], pmd_data)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.92092, 0.91012, 0.90059], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-deg2rad.([-30.0, -150.4, 89.8]), Inf) <= 0.1
        end
    end


    # @testset "2w transformer ac pf yy - banked transformers" begin
    #     file = "../test/data/opendss/ut_trans_2w_yy_bank.dss"
    #     pmd1 = PMD.parse_file(file)
    #     pmd2 = PMD.parse_file(file; bank_transformers=false)
    #     result1 = run_ac_mc_pf(pmd1, ipopt_solver)
    #     result2 = run_ac_mc_pf(pmd2, ipopt_solver)
    #
    #     @test result1["termination_status"] == PMs.LOCALLY_SOLVED
    #     @test result2["termination_status"] == PMs.LOCALLY_SOLVED
    #     @test result1["solution"]["bus"] == result2["solution"]["bus"]
    #     @test result1["solution"]["gen"] == result2["solution"]["gen"]
    #
    #     dss = PMD.parse_dss(file)
    #     PMD.parse_dss_with_dtypes!(dss, ["line", "load", "transformer"])
    #     trans = PMD._create_transformer(dss["transformer"][1]["name"]; PMD._to_sym_keys(dss["transformer"][1])...)
    #     @test all(trans["%rs"] .== [1.0, 2.0])
    # end

    @testset "three winding transformer pf" begin
        @testset "3w transformer ac pf dyy - all non-zero"  begin
            file = "../test/data/opendss/ut_trans_3w_dyy_1.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_mc_pf(pmd_data, ipopt_solver, multiconductor=true)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.9318, 0.88828, 0.88581], Inf) <= 1.5E-5
            @test norm(va(sol, pmd_data, "3")-[30.1, -90.7, 151.2], Inf) <= 0.1
        end

        @testset "3w transformer ac pf dyy - some non-zero" begin
            file = "../test/data/opendss/ut_trans_3w_dyy_2.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_mc_pf(pmd_data, ipopt_solver, multiconductor=true)
            #@test isapprox(vm(sol, pmd_data, "3"), [0.93876, 0.90227, 0.90454], atol=1E-4)
            @test norm(vm(sol, pmd_data, "3")-[0.93876, 0.90227, 0.90454], Inf) <= 1.5E-5
            @test norm(va(sol, pmd_data, "3")-[31.6, -88.8, 153.3], Inf) <= 0.1
        end

        @testset "3w transformer ac pf dyy - all zero" begin
            file = "../test/data/opendss/ut_trans_3w_dyy_3.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_mc_pf(pmd_data, ipopt_solver, multiconductor=true)
            solution_identify!(sol["solution"], pmd_data)
            @test norm(sol["solution"]["bus"]["3"]["vm"]-[0.97047, 0.93949, 0.946], Inf) <= 1.5E-5
            @test norm(sol["solution"]["bus"]["3"]["va"]-rad2deg.([30.6, -90.0, 151.9]), Inf) <= 0.1
        end

        @testset "3w transformer ac pf dyy - %loadloss=0" begin
            file = "../test/data/opendss/ut_trans_3w_dyy_3_loadloss.dss"
            pmd_data = PMD.parse_file(file)
            sol = PMD.run_ac_mc_pf(pmd_data, ipopt_solver, multiconductor=true)
            @test haskey(sol["solution"]["bus"], "10")
            @test norm(vm(sol, pmd_data, "3")-[0.969531, 0.938369, 0.944748], Inf) <= 1.5E-5
            @test norm(va(sol, pmd_data, "3")-[30.7, -90.0, 152.0], Inf) <= 0.1
        end
    end

    @testset "oltc tests" begin
        @testset "2w transformer acp opf_oltc yy" begin
            file = "../test/data/opendss/ut_trans_2w_yy_oltc.dss"
            pmd_data = PMD.parse_file(file)
            # free the taps
            pmd_data["transformer"]["1"]["fixed"] = zeros(Bool, 3)
            pmd_data["transformer"]["2"]["fixed"] = zeros(Bool, 3)
            pm = PMs.instantiate_model(pmd_data, PMs.ACPPowerModel, PMD.build_mc_opf_oltc, ref_extensions=[PMD.ref_add_arcs_trans!], multiconductor=true)
            sol = PMs.optimize_model!(pm, optimizer=ipopt_solver)
            # check that taps are set as to boost the voltage in the branches as much as possible;
            # this is trivially optimal if the voltage bounds are not binding
            # and without significant shunts (both branch and transformer)
            @test norm(tap(1,pm)-[0.95, 0.95, 0.95], Inf) <= 1E-4
            @test norm(tap(2,pm)-[1.05, 1.05, 1.05], Inf) <= 1E-4
            # then check whether voltage is what OpenDSS expects for those values
            @test norm(vm(sol, pmd_data, "3")-[1.0352, 1.022, 1.0142], Inf) <= 1E-4
            @test norm(va(sol, pmd_data, "3")-[-0.1, -120.4, 119.8], Inf) <= 0.1
        end
    end
end
