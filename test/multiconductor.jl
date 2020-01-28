# TESTLOG = Memento.getlogger(PowerModels)
calc_vm_acr(result, id) =   abs.(result["solution"]["bus"][id]["vr"] +im* result["solution"]["bus"][id]["vi"])
calc_va_acr(result, id) = angle.(result["solution"]["bus"][id]["vr"] +im* result["solution"]["bus"][id]["vi"])


# "an example of building a multi-phase model in an extention package"
# function build_tp_opf(pm::PMs.AbstractPowerModel)
#     PMD.variable_mc_voltage(pm)
#     PMD.variable_mc_branch_flow(pm)
#     PMD.variable_mc_transformer_flow(pm)
#     PMD.variable_mc_generation(pm)
#
#     PMD.constraint_mc_model_voltage(pm)
#
#     for i in PMs.ids(pm, :ref_buses)
#         PMD.constraint_mc_theta_ref(pm, i)
#     end
#
#     for i in PMs.ids(pm, :bus)
#         PMD.constraint_mc_power_balance(pm, i)
#     end
#
#     for i in PMs.ids(pm, :branch)
#         PMD.constraint_mc_ohms_yt_from(pm, i)
#         PMD.constraint_mc_ohms_yt_to(pm, i)
#
#         PMD.constraint_mc_voltage_angle_difference(pm, i)
#
#         PMD.constraint_mc_thermal_limit_from(pm, i)
#         PMD.constraint_mc_thermal_limit_to(pm, i)
#     end
#
#     for i in PMs.ids(pm, :transformer)
#         PMD.constraint_mc_trans(pm, i)
#     end
#
#     PMs.objective_min_fuel_cost(pm)
# end


@testset "test multi-conductor" begin

    @testset "json parser" begin
        mc_data = build_mc_data!("../test/data/pti/parser_test_defaults.raw")
        mc_json_string = PowerModels.parse_json(JSON.json(mc_data))
        @test mc_data == mc_json_string

        io = PipeBuffer()
        JSON.print(io, mc_data)
        mc_json_file = PowerModels.parse_file(io)
        @test mc_data == mc_json_file

        mc_strg_data = build_mc_data!("../test/data/matpower/case5_strg.m")
        mc_strg_json_string = PowerModels.parse_json(JSON.json(mc_strg_data))
        @test mc_strg_data == mc_strg_json_string

        # test that non-multiconductor json still parses, pti_json_file will result in error if fails
        pti_data = PowerModels.parse_file("../test/data/pti/parser_test_defaults.raw")
        io = PipeBuffer()
        JSON.print(io, pti_data)
        pti_json_file = PowerModels.parse_file(io)
        @test pti_data == pti_json_file

        mc_data = build_mc_data!("../test/data/matpower/case5.m")
        mc_data["gen"]["1"]["pmax"] = PMD.MultiConductorVector([Inf, Inf, Inf])
        mc_data["gen"]["1"]["qmin"] = PMD.MultiConductorVector([-Inf, -Inf, -Inf])
        mc_data["gen"]["1"]["bool_test"] = PMD.MultiConductorVector([true, true, false])
        mc_data["gen"]["1"]["string_test"] = PMD.MultiConductorVector(["a", "b", "c"])
        mc_data["branch"]["1"]["br_x"][1,2] = -Inf
        mc_data["branch"]["1"]["br_x"][1,3] = Inf

        mc_data_json = PowerModels.parse_json(JSON.json(mc_data))
        @test mc_data_json == mc_data

        mc_data["gen"]["1"]["nan_test"] = PMD.MultiConductorVector([0, NaN, 0])
        mc_data_json = PowerModels.parse_json(JSON.json(mc_data))
        @test isnan(mc_data_json["gen"]["1"]["nan_test"][2])
    end

    @testset "idempotent unit transformation" begin
        @testset "5-bus replicate case" begin
            mp_data = build_mc_data!("../test/data/matpower/case5_dc.m")

            PowerModels.make_mixed_units!(mp_data)
            PowerModels.make_per_unit!(mp_data)

            @test InfrastructureModels.compare_dict(mp_data, build_mc_data!("../test/data/matpower/case5_dc.m"))
        end
        @testset "24-bus replicate case" begin
            mp_data = build_mc_data!("../test/data/matpower/case24.m")

            PowerModels.make_mixed_units!(mp_data)
            PowerModels.make_per_unit!(mp_data)

            @test InfrastructureModels.compare_dict(mp_data, build_mc_data!("../test/data/matpower/case24.m"))
        end
    end


    @testset "topology processing" begin
        @testset "7-bus replicate status case" begin
            mp_data = build_mc_data!("../test/data/matpower/case7_tplgy.m")
            PowerModels.propagate_topology_status!(mp_data)
            PowerModels.select_largest_component!(mp_data)

            active_buses = Set(["4", "5", "7"])
            active_branches = Set(["8"])
            active_dclines = Set(["3"])

            for (i,bus) in mp_data["bus"]
                if i in active_buses
                    @test bus["bus_type"] != 4
                else
                    @test bus["bus_type"] == 4
                end
            end

            for (i,branch) in mp_data["branch"]
                if i in active_branches
                    @test branch["br_status"] == 1
                else
                    @test branch["br_status"] == 0
                end
            end

            for (i,dcline) in mp_data["dcline"]
                if i in active_dclines
                    @test dcline["br_status"] == 1
                else
                    @test dcline["br_status"] == 0
                end
            end
        end

    end


    @testset "test multi-conductor acp opf" begin
        @testset "3-bus 3-conductor case" begin
            mp_data = build_mc_data!("../test/data/matpower/case3.m", conductors=3)
            result = PMD.run_mc_opf(mp_data, PowerModels.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 47267.9; atol = 1e-1)

            for c in 1:mp_data["conductors"]
                @test isapprox(result["solution"]["gen"]["1"]["pg"][c], 1.58067; atol = 1e-3)
                @test isapprox(result["solution"]["bus"]["2"]["va"][c], 0.12669; atol = 1e-3)
            end
        end

        @testset "3-bus 3-conductor unbalanced case" begin
            # TESTS for add_setpoint! bug on multiconductor systems
            mp_data = build_mc_data!("../test/data/matpower/case3.m", conductors=3)
            for load in values(mp_data["load"])
                load["pd"][2] /= 2
                load["pd"][3] /= 3
            end
            result = PMD.run_mc_opf(mp_data, PowerModels.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 17826.8; atol = 1e-1)

            for (c, pg, va) in zip(1:mp_data["conductors"], [1.70667, 0.53344, 0.21976], [0.02996, 0.18645, 0.19194])
                @test isapprox(result["solution"]["gen"]["1"]["pg"][c], pg; atol = 1e-3)
                @test isapprox(result["solution"]["bus"]["2"]["va"][c], va; atol = 1e-3)
            end
        end

        @testset "3-bus 3-conductor case with theta_ref=pi" begin
            mp_data = build_mc_data!("../test/data/matpower/case3.m", conductors=3)
            pm = PowerModels.instantiate_model(mp_data, PowerModels.ACRPowerModel, PMD.build_mc_opf, multiconductor=true, ref_extensions=[PMD.ref_add_arcs_trans!])
            result = PowerModels.optimize_model!(pm, optimizer=ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 47267.9; atol = 1e-1)

            for c in 1:mp_data["conductors"]
                @test isapprox(result["solution"]["gen"]["1"]["pg"][c], 1.58067; atol = 1e-3)
                va = calc_va_acr(result, "2")
                @test isapprox(va, 0.12669; atol = 1e-3)
            end
        end

        @testset "5-bus 5-conductor case" begin
            mp_data = build_mc_data!("../test/data/matpower/case5.m", conductors=5)

            result = PMD.run_mc_opf(mp_data, PowerModels.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 91345.5; atol = 1e-1)
            for c in 1:mp_data["conductors"]
                @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.4; atol = 1e-3)
                @test isapprox(result["solution"]["bus"]["2"]["va"][c], -0.00103692; atol = 1e-5)
            end
        end

        @testset "30-bus 3-conductor case" begin
            mp_data = build_mc_data!("../test/data/matpower/case30.m", conductors=3)

            result = PMD.run_mc_opf(mp_data, PowerModels.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 614.905; atol = 1e-1)

            for c in 1:mp_data["conductors"]
                @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  2.18839; atol = 1e-3)
                @test isapprox(result["solution"]["bus"]["2"]["va"][c], -0.071759; atol = 1e-4)
            end
        end
    end


    @testset "test multi-conductor opf variants" begin
        mp_data = build_mc_data!("../test/data/matpower/case5_dc.m")

        @testset "ac 5-bus case" begin
            result = PMD.run_mc_opf(mp_data, PowerModels.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 54468.5; atol = 1e-1)
            for c in 1:mp_data["conductors"]
                @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.4; atol = 1e-3)
                @test isapprox(result["solution"]["bus"]["2"]["va"][c], -0.0139117; atol = 1e-4)
            end
        end


        # @testset "iv 5-bus case" begin
        #     result = PMD._run_mc_opf_iv(mp_data, PowerModels.IVRPowerModel, ipopt_solver)
        #
        #     @test result["termination_status"] == LOCALLY_SOLVED
        #     @test isapprox(result["objective"], 54468.5; atol = 1e-1)
        #     for c in 1:mp_data["conductors"]
        #         @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.4; atol = 1e-3)
        #         @test isapprox(result["solution"]["bus"]["2"]["va"][c], -0.0139117; atol = 1e-4)
        #     end
        # end

        @testset "dc 5-bus case" begin
            result = PMD.run_mc_opf(mp_data, PowerModels.DCPPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 54272.7; atol = 1e-1)
            for c in 1:mp_data["conductors"]
                @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.4; atol = 1e-3)
                @test isapprox(result["solution"]["bus"]["2"]["va"][c], -0.0135206; atol = 1e-4)
            end
        end

        @testset "soc 5-bus case" begin
            result = PMD.run_mc_opf(mp_data, PowerModels.SOCWRPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 46314.1; atol = 1e-1)
            for c in 1:mp_data["conductors"]
                vm = calc_vm_w(result, "2")
                @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.4; atol = 1e-3)
                @test isapprox(vm[c],  1.08578; atol = 1e-3)
            end
        end

    end


    @testset "test multi-conductor uc opf variants" begin
        mp_uc_data = build_mc_data!("../test/data/matpower/case5_uc.m")
        mp_strg_data = build_mc_data!("../test/data/matpower/case5_strg.m")

        @testset "ac 5-bus uc case" begin
            result = PMD._run_mc_ucopf(mp_uc_data, PowerModels.ACPPowerModel, juniper_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 54810.0; atol = 1e-1)
            @test isapprox(result["solution"]["gen"]["4"]["gen_status"], 0.0, atol=1e-6)
        end

        @testset "ac 5-bus storage case" begin
            result = PMD._run_mc_ucopf(mp_strg_data, PowerModels.ACPPowerModel, juniper_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 52633.8; atol = 1e-1)
        end


        @testset "dc 5-bus uc case" begin
            result = PMD._run_mc_ucopf(mp_uc_data, PowerModels.DCPPowerModel, cbc_solver)

            @test result["termination_status"] == OPTIMAL
            @test isapprox(result["objective"], 52839.6; atol = 1e-1)
            @test isapprox(result["solution"]["gen"]["4"]["gen_status"], 0.0)
        end

        @testset "dc 5-bus storage case" begin
            result = PMD._run_mc_ucopf(mp_strg_data, PowerModels.DCPPowerModel, cbc_solver)

            @test result["termination_status"] == OPTIMAL
            @test isapprox(result["objective"], 52081.3; atol = 1e-1)
        end

    end


    @testset "dual variable case" begin

        @testset "test dc polar opf" begin
            mp_data = build_mc_data!("../test/data/matpower/case5.m")

            result = PMD.run_mc_opf(mp_data, PowerModels.DCPPowerModel, ipopt_solver, setting = Dict("output" => Dict("duals" => true)))

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 52839.6; atol = 1e0)


            for (i, bus) in result["solution"]["bus"]
                @test haskey(bus, "lam_kcl_r")
                @test !haskey(bus, "lam_kcl_i")

                for c in 1:mp_data["conductors"]
                    @test bus["lam_kcl_r"][c] >= -4000 && bus["lam_kcl_r"][c] <= 0
                end
            end
            for (i, branch) in result["solution"]["branch"]
                @test haskey(branch, "mu_sm_fr")
                @test haskey(branch, "mu_sm_to")

                for c in 1:mp_data["conductors"]
                    @test branch["mu_sm_fr"][c] >= -1 && branch["mu_sm_fr"][c] <= 6000
                    @test isapprox(branch["mu_sm_to"][c], 0.0; atol = 1e-2)
                end
            end

        end
    end


    @testset "test solution feedback" begin
        mp_data = build_mc_data!("../test/data/matpower/case5_asym.m")

        result = PMD.run_mc_opf(mp_data, PowerModels.ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 52655.7; atol = 1e0)

        PowerModels.update_data!(mp_data, result["solution"])

        @test !InfrastructureModels.compare_dict(mp_data, build_mc_data!("../test/data/matpower/case5_asym.m"))
    end


    @testset "test errors and warnings" begin
        mp_data_2p = PowerModels.parse_file("../test/data/matpower/case3.m")
        mp_data_3p = PowerModels.parse_file("../test/data/matpower/case3.m")

        PMD.make_multiconductor!(mp_data_2p, 2)
        PMD.make_multiconductor!(mp_data_3p, 3)

        @test_throws(TESTLOG, ErrorException, PowerModels.update_data!(mp_data_2p, mp_data_3p))
        @test_throws(TESTLOG, ErrorException, PowerModels._check_keys(mp_data_3p, ["load"]))

        # check_cost_functions
        mp_data_3p["gen"]["1"]["model"] = 1
        mp_data_3p["gen"]["1"]["ncost"] = 1
        mp_data_3p["gen"]["1"]["cost"] = [0.0, 1.0, 0.0]
        @test_throws(TESTLOG, ErrorException, PowerModels.correct_cost_functions!(mp_data_3p))

        mp_data_3p["gen"]["1"]["cost"] = [0.0, 0.0]
        @test_throws(TESTLOG, ErrorException, PowerModels.correct_cost_functions!(mp_data_3p))

        mp_data_3p["gen"]["1"]["ncost"] = 2
        mp_data_3p["gen"]["1"]["cost"] = [0.0, 1.0, 0.0, 2.0]
        @test_throws(TESTLOG, ErrorException, PowerModels.correct_cost_functions!(mp_data_3p))

        mp_data_3p["gen"]["1"]["model"] = 2
        @test_throws(TESTLOG, ErrorException, PowerModels.correct_cost_functions!(mp_data_3p))

        Memento.setlevel!(TESTLOG, "info")

        mp_data_3p["gen"]["1"]["model"] = 3
        @test_warn(TESTLOG, "Skipping cost model of type 3 in per unit transformation", PowerModels.make_mixed_units!(mp_data_3p))
        @test_warn(TESTLOG, "Skipping cost model of type 3 in per unit transformation", PowerModels.make_per_unit!(mp_data_3p))
        @test_warn(TESTLOG, "Unknown cost model of type 3 on generator 1", PowerModels.correct_cost_functions!(mp_data_3p))

        mp_data_3p["gen"]["1"]["model"] = 1
        mp_data_3p["gen"]["1"]["cost"][3] = 3000
        @test_warn(TESTLOG, "pwl x value 3000.0 is outside the bounds 0.0-60.0 on generator 1", PowerModels.correct_cost_functions!(mp_data_3p))

        @test_nowarn PowerModels.correct_voltage_angle_differences!(mp_data_3p)

        mp_data_2p["branch"]["1"]["angmin"] = [-pi, 0]
        mp_data_2p["branch"]["1"]["angmax"] = [ pi, 0]

        @test_warn(TESTLOG, "this code only supports angmin values in -90 deg. to 90 deg., tightening the value on branch 1, conductor 1 from -180.0 to -60.0 deg.",
            PowerModels.correct_voltage_angle_differences!(mp_data_2p))

        mp_data_2p["branch"]["1"]["angmin"] = [-pi, 0]
        mp_data_2p["branch"]["1"]["angmax"] = [ pi, 0]

        @test_warn(TESTLOG, "angmin and angmax values are 0, widening these values on branch 1, conductor 2 to +/- 60.0 deg.",
            PowerModels.correct_voltage_angle_differences!(mp_data_2p))

        mp_data_2p["branch"]["1"]["angmin"] = [-pi, 0]
        mp_data_2p["branch"]["1"]["angmax"] = [ pi, 0]

        @test_warn(TESTLOG, "this code only supports angmax values in -90 deg. to 90 deg., tightening the value on branch 1, conductor 1 from 180.0 to 60.0 deg.",
            PowerModels.correct_voltage_angle_differences!(mp_data_2p))

        @test_warn(TESTLOG, "skipping network that is already multiconductor", PMD.make_multiconductor!(mp_data_3p, 3))

        mp_data_3p["load"]["1"]["pd"] = mp_data_3p["load"]["1"]["qd"] = [0, 0, 0]
        mp_data_3p["shunt"]["1"] = Dict("gs"=>[0,0,0], "bs"=>[0,0,0], "status"=>1, "shunt_bus"=>1, "index"=>1)

        Memento.TestUtils.@test_log(TESTLOG, "info", "deactivating load 1 due to zero pd and qd", PowerModels.propagate_topology_status!(mp_data_3p))
        @test mp_data_3p["load"]["1"]["status"] == 0
        @test mp_data_3p["shunt"]["1"]["status"] == 0

        mp_data_3p["dcline"]["1"]["loss0"][2] = -1.0
        @test_warn(TESTLOG, "this code only supports positive loss0 values, changing the value on dcline 1, conductor 2 from -100.0 to 0.0", PowerModels.correct_dcline_limits!(mp_data_3p))

        mp_data_3p["dcline"]["1"]["loss1"][2] = -1.0
        @test_warn(TESTLOG, "this code only supports positive loss1 values, changing the value on dcline 1, conductor 2 from -1.0 to 0.0", PowerModels.correct_dcline_limits!(mp_data_3p))

        @test mp_data_3p["dcline"]["1"]["loss0"][2] == 0.0
        @test mp_data_3p["dcline"]["1"]["loss1"][2] == 0.0

        mp_data_3p["dcline"]["1"]["loss1"][2] = 100.0
        @test_warn(TESTLOG, "this code only supports loss1 values < 1, changing the value on dcline 1, conductor 2 from 100.0 to 0.0", PowerModels.correct_dcline_limits!(mp_data_3p))

        delete!(mp_data_3p["branch"]["1"], "tap")
        @test_warn(TESTLOG, "branch found without tap value, setting a tap to 1.0", PowerModels.correct_transformer_parameters!(mp_data_3p))

        delete!(mp_data_3p["branch"]["1"], "shift")
        @test_warn(TESTLOG, "branch found without shift value, setting a shift to 0.0", PowerModels.correct_transformer_parameters!(mp_data_3p))

        mp_data_3p["branch"]["1"]["tap"][2] = -1.0
        @test_warn(TESTLOG, "branch found with non-positive tap value of -1.0, setting a tap to 1.0", PowerModels.correct_transformer_parameters!(mp_data_3p))

        mp_data_3p["branch"]["1"]["rate_a"][2] = -1.0
        @test_throws(TESTLOG, ErrorException, PowerModels.correct_thermal_limits!(mp_data_3p))
        mp_data_3p["branch"]["1"]["rate_a"][2] = 1.0

        mp_data_3p["branch"]["4"] = deepcopy(mp_data_3p["branch"]["1"])
        mp_data_3p["branch"]["4"]["f_bus"] = mp_data_3p["branch"]["1"]["t_bus"]
        mp_data_3p["branch"]["4"]["t_bus"] = mp_data_3p["branch"]["1"]["f_bus"]
        @test_warn(TESTLOG, "reversing the orientation of branch 1 (1, 3) to be consistent with other parallel branches", PowerModels.correct_branch_directions!(mp_data_3p))
        @test mp_data_3p["branch"]["4"]["f_bus"] == mp_data_3p["branch"]["1"]["f_bus"]
        @test mp_data_3p["branch"]["4"]["t_bus"] == mp_data_3p["branch"]["1"]["t_bus"]

        mp_data_3p["gen"]["1"]["vg"][2] = 2.0
        @test_warn(TESTLOG, "the conductor 2 voltage setpoint on generator 1 does not match the value at bus 1", PowerModels.check_voltage_setpoints(mp_data_3p))

        mp_data_3p["dcline"]["1"]["vf"][2] = 2.0
        @test_warn(TESTLOG, "the conductor 2 from bus voltage setpoint on dc line 1 does not match the value at bus 1", PowerModels.check_voltage_setpoints(mp_data_3p))

        mp_data_3p["dcline"]["1"]["vt"][2] = 2.0
        @test_warn(TESTLOG, "the conductor 2 to bus voltage setpoint on dc line 1 does not match the value at bus 2", PowerModels.check_voltage_setpoints(mp_data_3p))


        Memento.setlevel!(TESTLOG, "error")

        @test_throws(TESTLOG, ErrorException, PowerModels.run_ac_opf(mp_data_3p, ipopt_solver))

        mp_data_3p["branch"]["1"]["f_bus"] = mp_data_3p["branch"]["1"]["t_bus"] = 1
        @test_throws(TESTLOG, ErrorException, PowerModels.check_branch_loops(mp_data_3p))

        mp_data_3p["conductors"] = 0
        @test_throws(TESTLOG, ErrorException, PowerModels.check_conductors(mp_data_3p))
    end

    @testset "multiconductor extensions" begin
        mp_data = build_mc_data!("../test/data/matpower/case3.m")
        pm = PowerModels.instantiate_model(mp_data, PowerModels.ACPPowerModel, build_mc_opf; multiconductor=true, ref_extensions=[PMD.ref_add_arcs_trans!])

        @test length(PMs.var(pm, pm.cnw)) == 13

        @test haskey(PMs.var(pm, pm.cnw), :vm)
        @test length(PMs.var(pm, pm.cnw, :vm)) == 3

        @test PowerModels.ref(pm, pm.cnw, :bus, 1)["bus_i"] == 1
        @test PowerModels.ref(pm, :bus, 1)["vmax"][1] == 1.1

        @test PowerModels.ismulticonductor(pm)
        @test PowerModels.ismulticonductor(pm, pm.cnw)

        @test length(PowerModels.nw_ids(pm)) == 1
    end

    @testset "multiconductor operations" begin
        mp_data = build_mc_data!("../test/data/matpower/case3.m")

        a, b, c, d = mp_data["branch"]["1"]["br_r"], mp_data["branch"]["1"]["br_x"], mp_data["branch"]["1"]["b_fr"], mp_data["branch"]["1"]["b_to"]
        c = diag(c)
        d = diag(d)
        e = PMD.MultiConductorVector([0.225, 0.225, 0.225, 0.225])
        angs_rad = mp_data["branch"]["1"]["angmin"]

        # Transpose
        @test all(a' .== a)
        @test all(c' .== [0.225, 0.225, 0.225]')

        # Basic Math (Matrices)
        x = a + b
        y = a - b
        z = a * b
        w = a / b

        @test all(x.values - [0.685 0.0 0.0; 0.0 0.685 0.0; 0.0 0.0 0.685] .<= 1e-12)
        @test all(y.values - [-0.555 0.0 0.0; 0.0 -0.555 0.0; 0.0 0.0 -0.555] .<= 1e-12)
        @test all(z.values - [0.0403 0.0 0.0; 0.0 0.0403 0.0; 0.0 0.0 0.0403] .<= 1e-12)
        @test all(w.values - [0.104839 0.0 0.0; 0.0 0.104839 0.0; 0.0 0.0 0.104839] .<= 1e-12)

        @test isa(x, PMD.MultiConductorMatrix)
        @test isa(y, PMD.MultiConductorMatrix)
        @test isa(z, PMD.MultiConductorMatrix)
        @test isa(w, PMD.MultiConductorMatrix)

        # Basic Math Vectors
        x = c + d
        y = c - d
        z = c^2
        w = 1 ./ c
        u = 1 .* d

        @test all(x.values - [0.45, 0.45, 0.45] .<= 1e-12)
        @test all(y.values - [0.0, 0.0, 0.0] .<= 1e-12)
        @test all(c .* d - [0.050625, 0.050625, 0.050625] .<= 1e-12)
        @test all(c ./ d - [1.0, 1.0, 1.0] .<= 1e-12)
        @test all(z.values - [0.050625, 0.050625, 0.050625] .<= 1e-12)
        @test all(w.values - [4.444444444444445, 4.444444444444445, 4.444444444444445] .<= 1e-12)
        @test all(u.values - d.values .<= 1e-12)

        @test isa(x, PMD.MultiConductorVector)
        @test isa(y, PMD.MultiConductorVector)
        @test isa(z, PMD.MultiConductorVector)
        @test isa(w, PMD.MultiConductorVector)
        @test isa(u, PMD.MultiConductorVector)

        # Broadcasting
        @test all(a .+ c - [0.29   0.225  0.225; 0.225  0.29   0.225; 0.225  0.225  0.29] .<= 1e-12)
        @test all(c .+ b - [0.845  0.225  0.225; 0.225  0.845  0.225; 0.225  0.225  0.845] .<= 1e-12)
        @test all(a.values .+ c - [0.29   0.225  0.225; 0.225  0.29   0.225; 0.225  0.225  0.29] .<= 1e-12)
        @test all(c .+ b.values - [0.845  0.225  0.225; 0.225  0.845  0.225; 0.225  0.225  0.845] .<= 1e-12)

        # Custom Functions
        @test PMD.conductors(c) == 3
        @test PMD.conductors(a) == 3
        @test all(size(a) == (3,3))
        @test isa(JSON.lower(a), Dict)
        @test all(JSON.lower(a)["values"] == a.values)
        @test !isapprox(d, e)
        @test PMD.conductor_value(a, 1, 1) == a[1,1]

        # diagm
        @test all(LinearAlgebra.diagm(0 => c).values .== [0.225 0.0 0.0; 0.0 0.225 0.0; 0.0 0.0 0.225])

        # rad2deg/deg2rad
        angs_deg = rad2deg(angs_rad)
        angs_deg_rad = deg2rad(angs_deg)
        @test all(angs_deg.values - [-30.0, -30.0, -30.0] .<= 1e-12)
        @test all(angs_deg_rad.values - angs_rad.values .<= 1e-12)

        @test isa(angs_deg, PMD.MultiConductorVector)
        @test isa(deg2rad(angs_deg), PMD.MultiConductorVector)

        a_rad = rad2deg(a)
        @test all(a_rad.values - [3.72423 0.0 0.0; 0.0 3.72423 0.0; 0.0 0.0 3.72423] .<= 1e-12)
        @test isa(rad2deg(a), PMD.MultiConductorMatrix)
        @test isa(deg2rad(a), PMD.MultiConductorMatrix)

        Memento.setlevel!(TESTLOG, "warn")
        @test_nowarn show(devnull, a)

        @test_nowarn a[1, 1] = 9.0
        @test a[1,1] == 9.0

        @test_nowarn PowerModels.summary(devnull, mp_data)

        Memento.setlevel!(TESTLOG, "error")

        # Test broadcasting edge-case
        v = ones(Real, 3)
        mcv = PMD.MultiConductorVector(v)
        @test all(floor.(mcv) .+ mcv .== PMD.MultiConductorVector(floor.(v) .+ v))

        m = LinearAlgebra.diagm(0 => v)
        mcm = PMD.MultiConductorMatrix(m)
        @test all(floor.(mcm) .+ mcm .== PMD.MultiConductorMatrix(floor.(m) .+ m))
    end
end
