@info "running optimal power flow (opf) tests"

@testset "test opf" begin
    @testset "test matpower opf" begin
        case5 = PM.parse_file("../test/data/matpower/case5.m")
        case5_strg = PM.parse_file("$(pms_path)/test/data/matpower/case5_strg.m")
        case30 = PM.parse_file("../test/data/matpower/case30.m")

        make_multiconductor!(case5, 3)
        make_multiconductor!(case5_strg, 3)
        make_multiconductor!(case30, 3)

        @testset "5-bus matpower acp opf" begin
            result = solve_mc_opf(case5, ACPUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 45522.096; atol=1e-1)

            @test all(isapprox.(result["solution"]["gen"]["1"]["pg"],  0.3999999; atol=1e-3))
            @test all(isapprox.(result["solution"]["bus"]["2"]["va"], -0.0538204; atol=1e-5))
        end

        @testset "5-bus matpower acr opf" begin
            result = solve_mc_opf(case5, ACRUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 45522.096; atol=1e-1)

            calc_va(id) = atan.(result["solution"]["bus"][id]["vi"], result["solution"]["bus"][id]["vr"])
            @test all(isapprox.(result["solution"]["gen"]["1"]["pg"],  0.3999999; atol=1e-3))
            @test all(isapprox.(calc_va("2"), -0.0538204; atol=1e-5))
        end

        @testset "5-bus matpower mn acp mld" begin
            case5_mn = InfrastructureModels.replicate(case5, 3, PowerModelsDistribution._pmd_math_global_keys)
            result = solve_mn_mc_opf(case5_mn, ACPUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 45522.096*3; atol=1e-1)

            @test all(isapprox.(result["solution"]["nw"]["1"]["gen"]["1"]["pg"],  0.3999999; atol=1e-3))
            @test all(isapprox.(result["solution"]["nw"]["3"]["bus"]["2"]["va"], -0.0538204; atol=1e-5))
        end

        @testset "5-bus matpower mn acr mld" begin
            case5_mn = InfrastructureModels.replicate(case5, 3, PowerModelsDistribution._pmd_math_global_keys)
            result = solve_mn_mc_opf(case5_mn, ACRUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 45522.096*3; atol=1e-1)

            calc_va(id) = atan.(result["solution"]["nw"]["2"]["bus"][id]["vi"], result["solution"]["nw"]["2"]["bus"][id]["vr"])
            @test all(isapprox.(result["solution"]["nw"]["3"]["gen"]["1"]["pg"], 0.3999999; atol=1e-3))
            @test all(isapprox.(calc_va("2"), -0.0538204; atol=1e-5))
        end

        @testset "5-bus storage matpower mn acr mld" begin
            case5_strg_mn = InfrastructureModels.replicate(case5_strg, 3, PowerModelsDistribution._pmd_math_global_keys)
            result = solve_mn_mc_opf(case5_strg_mn, ACRUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test result["objective"] >= 45522.096*3
        end

        @testset "30-bus matpower acp opf" begin
            result = solve_mc_opf(case30, ACPUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 614.007; atol=1e-1)

            @test all(isapprox.(result["solution"]["gen"]["1"]["pg"],  2.192189; atol=1e-3))
            @test all(isapprox.(result["solution"]["bus"]["2"]["va"], -0.071853; atol=1e-4))
        end

        @testset "30-bus matpower acr opf" begin
            result = solve_mc_opf(case30, ACRUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 614.007; atol=1e-1)

            calc_va(id) = atan.(result["solution"]["bus"][id]["vi"], result["solution"]["bus"][id]["vr"])
            @test all(isapprox.(result["solution"]["gen"]["1"]["pg"],  2.192189; atol=1e-3))
            @test all(isapprox.(calc_va("2"), -0.071853; atol=1e-4))
        end

        @testset "30-bus matpower dcp opf" begin
            result = solve_mc_opf(case30, DCPUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 566.112; atol=1e-1)
        end

        @testset "30-bus matpower nfa opf" begin
            result = solve_mc_opf(case30, NFAUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 458.006; atol=1e-1)
        end
    end

    @testset "test dropped phases opf" begin
        case4_phase_drop = parse_file("../test/data/opendss/case4_phase_drop.dss")
        case5_phase_drop = parse_file("../test/data/opendss/case5_phase_drop.dss")

        @testset "4-bus phase drop acp opf" begin
            result = solve_mc_opf(case4_phase_drop, ACPUPowerModel, ipopt_solver, make_si=false)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0182595; atol=1e-4)

            @test all(isapprox.(result["solution"]["voltage_source"]["source"]["pg"], [5.06513e-5, 6.0865e-5, 7.1119e-5]; atol=1e-7))
            @test isapprox(result["solution"]["bus"]["loadbus1"]["vm"][1], 0.98995; atol=1.5e-4)
        end

        @testset "4-bus phase drop acr opf" begin
            result = solve_mc_opf(case4_phase_drop, ACRUPowerModel, ipopt_solver; make_si=false)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0182595; atol=1e-4)

            @test all(isapprox.(result["solution"]["voltage_source"]["source"]["pg"], [5.06513e-5, 6.0865e-5, 7.1119e-5]; atol=1e-7))
            @test isapprox(calc_vm_acr(result, "loadbus1")[1], 0.98995; atol=1.5e-4)
        end

        @testset "5-bus phase drop acp opf" begin
            result = solve_mc_opf(case5_phase_drop, ACPUPowerModel, ipopt_solver; make_si=false)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0599389; atol=1e-4)

            @test all(isapprox.(result["solution"]["voltage_source"]["source"]["pg"], [0.000152, 0.000198, 0.000248]; atol=1e-6))
            @test all(isapprox.(result["solution"]["bus"]["midbus"]["vm"], [0.97351, 0.96490, 0.95646]; atol=1e-4))
        end

        @testset "5-bus phase drop acr opf" begin
            result = solve_mc_opf(case5_phase_drop, ACRUPowerModel, ipopt_solver; make_si=false)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.0599400; atol = 1e-4)

            @test all(isapprox.(result["solution"]["voltage_source"]["source"]["pg"], [0.00015236280779412599, 0.00019836795302238667, 0.0002486688793741932]; atol=1e-7))
            @test all(isapprox.(calc_vm_acr(result, "midbus"), [0.9735188343958152, 0.9649003198689144, 0.9564593296045091]; atol=1e-4))
        end

        @testset "5-bus phase drop dcp opf" begin
            result = solve_mc_opf(case5_phase_drop, DCPUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.054; atol=1e-4)
        end

        @testset "5-bus phase drop nfa opf" begin
            result = solve_mc_opf(case5_phase_drop, NFAUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.054; atol=1e-4)
        end

        @testset "test warm start from dss voltages export" begin
            eng = deepcopy(case5_phase_drop)
            add_voltage_starts!(eng, "../test/data/opendss/case5_voltage.csv")
            apply_voltage_bounds!(eng)

            for (_,bus) in eng["bus"]
                nt = length(bus["terminals"])
                @test haskey(bus, "vm_start") && length(bus["vm_start"]) == nt
                @test haskey(bus, "va_start") && length(bus["va_start"]) == nt
            end

            math = transform_data_model(eng)

            for (_,bus) in math["bus"]
                if startswith(bus["source_id"], "bus")
                    @test bus["vm_start"] <= bus["vmax"] && bus["vm_start"] >= bus["vmin"]
                end
            end

            pm = instantiate_mc_model(math, ACPUPowerModel, build_mc_opf)
            for (i,vm) in var(pm, nw_id_default, :vm)
                if startswith(math["bus"]["$i"]["source_id"], "bus")
                    @test all(JuMP.start_value(vm[t]) == math["bus"]["$i"]["vm_start"][idx] for (idx,t) in enumerate(math["bus"]["$i"]["terminals"]))
                end
            end

            for (i,va) in var(pm, nw_id_default, :va)
                if startswith(math["bus"]["$i"]["source_id"], "bus")
                    @test all(JuMP.start_value(va[t]) == math["bus"]["$i"]["va_start"][idx] for (idx,t) in enumerate(math["bus"]["$i"]["terminals"]))
                end
            end

            result = solve_mc_opf(eng, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)

            @test all(isapprox.(result["solution"]["bus"]["midbus"]["vm"], [0.97352, 0.9649, 0.95646]; atol=1e-4))
            @test all(isapprox.(result["solution"]["bus"]["midbus"]["va"], [-1.3, -121.3, 118.2]; atol=1e-1))
        end
    end

    @testset "test opendss opf" begin
        @testset "2-bus diagonal acp opf" begin
            pmd = parse_file("../test/data/opendss/case2_diag.dss")
            sol = solve_mc_opf(pmd, ACPUPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED

            @test all(isapprox.(sol["solution"]["bus"]["primary"]["vm"], 0.984377; atol=1e-4))
            @test all(isapprox.(sol["solution"]["bus"]["primary"]["va"], [0, -120, 120] .- 0.79; atol=0.2))

            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA), 0.0181409; atol=1e-5)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA), 0.0; atol=1e-4)
        end

        @testset "3-bus balanced acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_balanced.dss")
            sol = solve_mc_opf(pmd, ACPUPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED

            for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"], [0.0, -0.03, -0.07], [0.9959, 0.986973, 0.976605])
                @test all(isapprox.(sol["solution"]["bus"][bus]["va"], [0, -120, 120] .+ va; atol=0.01))
                @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=1e-4))
            end

            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA), 0.018276; atol=1e-6)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA), 0.008922; atol=1.2e-5)
        end

        @testset "3-bus unbalanced acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_unbalanced.dss")
            sol = solve_mc_opf(pmd, ACPUPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED

            for (bus, va, vm) in zip(["sourcebus", "primary", "loadbus"],
                                    [0.0, [-0.22, -0.11, 0.12], [-0.48, -0.24, 0.27]],
                                    [0.9959, [0.980937, 0.98936, 0.987039], [0.963546, 0.981757, 0.976779]])
                @test all(isapprox.(sol["solution"]["bus"][bus]["va"], [0, -120, 120] .+ va; atol=0.01))
                @test all(isapprox.(sol["solution"]["bus"][bus]["vm"], vm; atol=1e-5))
            end

            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA), 0.0214812; atol=1e-6)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA), 0.00927263; atol=1e-5)
        end

        @testset "3-bus unbalanced isc acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_balanced_isc.dss")
            sol = solve_mc_opf(pmd, ACPUPowerModel, ipopt_solver)

            @test sol["termination_status"] == LOCALLY_SOLVED
            @test isapprox(sol["objective"], 0.0185; atol=1e-4)
        end

        @testset "3-bus balanced pv acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_balanced_pv.dss")

            sol = solve_mc_opf(pmd, ACPUPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED
            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA) < 0.0
            @test sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA) < 0.005
            @test isapprox(sum(sol["solution"]["solar"]["pv1"]["pg"] * baseMVA), 0.0183685; atol=1e-4)
            @test isapprox(sum(sol["solution"]["solar"]["pv1"]["qg"] * baseMVA), 0.0048248; atol=1e-4)
        end

        @testset "3-bus unbalanced single-phase pv acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_unbalanced_1phase-pv.dss")
            sol = solve_mc_opf(pmd, ACPUPowerModel, ipopt_solver)

            @test sol["termination_status"] == LOCALLY_SOLVED

            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"]), 18.38728; atol=1e-3)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"]),  7.28903; atol=1e-3)

            @test all(isapprox.(sol["solution"]["solar"]["pv1"]["pg"], 1.9947; atol=1e-3))
            @test all(isapprox.(sol["solution"]["solar"]["pv1"]["qg"], 1.9259; atol=1e-3))
        end

        @testset "3-bus balanced capacitor acp opf" begin
            pmd = parse_file("../test/data/opendss/case3_balanced_cap.dss")
            sol = solve_mc_pf(pmd, ACPUPowerModel, ipopt_solver; make_si=false)

            @test sol["termination_status"] == LOCALLY_SOLVED

            @test all(abs.(sol["solution"]["bus"]["loadbus"]["vm"].-0.98588).<=1E-4)
            @test all(abs.(sol["solution"]["bus"]["primary"]["vm"].-0.99127).<=1E-4)
        end

        @testset "3w transformer nfa opf" begin
            mp_data = parse_file("../test/data/opendss/ut_trans_3w_dyy_1.dss")
            result = solve_mc_opf(mp_data, NFAUPowerModel, ipopt_solver)

            @test result["termination_status"] == LOCALLY_SOLVED
            @test isapprox(result["objective"], 0.616; atol=1e-3)

        end

        @testset "vm and va start with ivr and acr opf" begin

            pmd_eng = parse_file("../test/data/opendss/case3_unbalanced.dss")
            pmd_math = transform_data_model(pmd_eng)

            sol_ivr = solve_mc_opf(pmd_math, IVRUPowerModel, ipopt_solver)
            sol_acr = solve_mc_opf(pmd_math, ACRUPowerModel, ipopt_solver)

            pmd_math["bus"]["4"]["vm_start"] = pmd_math["bus"]["4"]["vm"]
            pmd_math["bus"]["4"]["va_start"] = pmd_math["bus"]["4"]["va"]
            pmd_math["bus"]["2"]["vm_start"] = [0.9959, 0.9959, 0.9959]
            pmd_math["bus"]["2"]["va_start"] = [0.00, -2.0944, 2.0944]

            sol_ivr_with_start = solve_mc_opf(pmd_math, IVRUPowerModel, ipopt_solver)
            sol_acr_with_start = solve_mc_opf(pmd_math, ACRUPowerModel, ipopt_solver)

            @test isapprox(sol_ivr["objective"], sol_ivr_with_start["objective"]; atol=1e-5)
            @test isapprox(sol_acr["objective"], sol_acr_with_start["objective"]; atol=1e-5)
        end

        @testset "assign start value per connection iteration" begin
            data = parse_file("../test/data/opendss/case3_unbalanced.dss"; data_model = MATHEMATICAL)
            #add a single-phase generator and assign a single-phase start value to it
            data["gen"]["2"] = Dict{String, Any}(
                                "pg_start"      => [-0.012],
                                "qg_start"      => [-0.006],
                                "model"         => 2,
                                "connections"   => [2],
                                "shutdown"      => 0.0,
                                "startup"       => 0.0,
                                "configuration" => WYE,
                                "name"          => "single_ph_generator",
                                "gen_bus"       => 3,
                                "pmax"          => [Inf],
                                "vbase"         => 0.23094,
                                "index"         => 2,
                                "cost"          => [0.5, 0.0],
                                "gen_status"    => 1,
                                "qmax"          => [Inf],
                                "qmin"          => [-Inf],
                                "pmin"          => [-Inf],
                                "ncost"         => 2
                                )
            sol = solve_mc_opf(data, ACPUPowerModel, ipopt_solver)
            @test sol["termination_status"] == LOCALLY_SOLVED
        end
        @testset "3-bus unbalanced fotp opf with yy transformer" begin
            pmd = parse_file("../test/data/opendss/ut_trans_2w_yy.dss")
            sol = solve_mc_opf(pmd, FOTPUPowerModel, ipopt_solver; make_si=false)
            @test sol["termination_status"] == LOCALLY_SOLVED
            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA), 0.467547; atol=4e-2)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA), 0.484327; atol=3e-2)
            @test all(isapprox.(sol["solution"]["bus"]["3"]["vm"], [0.87451, 0.8613, 0.85348]; atol=2e-2))
            @test all(isapprox.(sol["solution"]["bus"]["3"]["va"], [-0.1, -120.4, 119.8]; atol=2e-1))
        end
        @testset "3-bus unbalanced fotp opf with dy transformer" begin
            pmd = parse_file("../test/data/opendss/ut_trans_2w_dy_lag.dss")
            sol = solve_mc_opf(pmd, FOTPUPowerModel, ipopt_solver; make_si=false)
            @test sol["termination_status"] == LOCALLY_SOLVED
            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA), 0.467699; atol=2e-1)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA), 0.485553; atol=6e-2)
            @test all(isapprox.(sol["solution"]["bus"]["3"]["vm"], [0.92092, 0.91012, 0.90059]; atol=3e-1))
        end
        @testset "3-bus unbalanced fotp opf with voltage-dependent loads" begin
            pmd = parse_file("../test/data/opendss/case3_unbalanced_delta_loads.dss")
            sol = solve_mc_opf(pmd, FOTPUPowerModel, ipopt_solver; make_si=false)
            @test sol["termination_status"] == LOCALLY_SOLVED
            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA), 0.0420464; atol=1e-3)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA), 0.0181928; atol=1e-3)
            @test all(isapprox.(sol["solution"]["bus"]["loadbus"]["vm"], [0.94105, 0.95942, 0.95876]; atol=2e-3))
            @test all(isapprox.(sol["solution"]["bus"]["loadbus"]["va"], [-0.9, -120.3, 120.2]; atol=2e-1))
        end
        @testset "3-bus unbalanced fotr opf with yy transformer" begin
            pmd = parse_file("../test/data/opendss/ut_trans_2w_yy.dss")
            sol = solve_mc_opf(pmd, FOTRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test sol["termination_status"] == LOCALLY_SOLVED
            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA), 0.467547; atol=7e-2)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA), 0.484327; atol=9e-2)
            @test all(isapprox.(sol["solution"]["bus"]["3"]["vm"], [0.87451, 0.8613, 0.85348]; atol=2e-2))
            @test all(isapprox.(sol["solution"]["bus"]["3"]["va"], [-0.1, -120.4, 119.8]; atol=2e-1))
        end
        @testset "3-bus unbalanced fotr opf with dy transformer" begin
            pmd = parse_file("../test/data/opendss/ut_trans_2w_dy_lag.dss")
            sol = solve_mc_opf(pmd, FOTRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test sol["termination_status"] == LOCALLY_SOLVED
            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA), 0.467699; atol=2e-1)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA), 0.485553; atol=6e-2)
            @test all(isapprox.(sol["solution"]["bus"]["3"]["vm"], [0.92092, 0.91012, 0.90059]; atol=3e-1))
        end
        @testset "3-bus unbalanced fotr opf with voltage-dependent loads" begin
            pmd = parse_file("../test/data/opendss/case3_unbalanced_delta_loads.dss")
            sol = solve_mc_opf(pmd, FOTRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!], make_si=false)
            @test sol["termination_status"] == LOCALLY_SOLVED
            baseMVA = sol["solution"]["settings"]["sbase"] / sol["solution"]["settings"]["power_scale_factor"]
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["pg"] * baseMVA), 0.0420464; atol=1e-3)
            @test isapprox(sum(sol["solution"]["voltage_source"]["source"]["qg"] * baseMVA), 0.0181928; atol=1e-3)
            @test all(isapprox.(sol["solution"]["bus"]["loadbus"]["vm"], [0.94105, 0.95942, 0.95876]; atol=2e-3))
            @test all(isapprox.(sol["solution"]["bus"]["loadbus"]["va"], [-0.9, -120.3, 120.2]; atol=2e-1))
        end
    end
end
