@info "running generator configuration tests"

@testset "test generator configuration" begin
    # This test checks the generators are connected properly by comparing them
    # to equivalent constant-power loads. This is achieved by fixing their bounds.
    @testset "ACP/ACR tests" begin
        eng_1 = parse_file("../test/data/opendss/case3_delta_gens.dss")

        for (_,load) in eng_1["load"]
            load["model"] = POWER
        end

        for (_,line) in eng_1["line"]
            line["cm_ub"] = fill(1e4, size(line["cm_ub"])...)
        end

        eng_2 = deepcopy(eng_1)
        eng_2["load"] = Dict{String,Any}()
        eng_2["generator"] = Dict{String,Any}()
        for (id,load) in eng_1["load"]
            gen = Dict{String,Any}(
                "source_id" => load["source_id"],
                "configuration" => load["configuration"],
                "bus" => load["bus"],
                "connections" => load["connections"],
                "cost_pg_parameters" => [0, 0, 0],
                "control_mode" => FREQUENCYDROOP,
                "pg_lb" => -load["pd_nom"],
                "pg_ub" => -load["pd_nom"],
                "qg_lb" => -load["qd_nom"],
                "qg_ub" => -load["qd_nom"],
                "status" => ENABLED,
            )

            eng_2["generator"][id] = gen
        end

        # check ACP and ACR
        for form in [ACPUPowerModel, ACRUPowerModel, IVRUPowerModel]
            sol_1 = solve_mc_opf(eng_1, form, ipopt_solver)
            @test sol_1["termination_status"] == LOCALLY_SOLVED

            sol_2 = solve_mc_opf(eng_2, form, ipopt_solver)
            @test sol_2["termination_status"] == LOCALLY_SOLVED

            # check that gens are equivalent to the loads
            for (id,_) in eng_1["load"]
                pd_bus = sol_1["solution"]["load"][id]["pd"]
                qd_bus = sol_1["solution"]["load"][id]["qd"]
                pg_bus = sol_2["solution"]["generator"][id]["pg"]
                qg_bus = sol_2["solution"]["generator"][id]["qg"]
                @test isapprox(pd_bus, -pg_bus, atol=1E-5)
                @test isapprox(qd_bus, -qg_bus, atol=1E-5)
            end
        end
    end
    # IVR cannot be tested this way; it does not have explicit power variables,
    # so setting the upper and lower bound for the power to the same variable
    # will be enforced indirectly, which is not very stable.
    # In the end it worked; keep this as a backup in case Ipopt has convergence
    # issue in the future on IVR in the previous test.
    # @testset "IVR" begin
    #     pmd = PMD.parse_file("../test/data/opendss/case3_delta_gens.dss")
    #
    #     for (i,(id, load)) in enumerate(pmd["load"])
    #         load = pmd_1["load"][id]
    #         gen = deepcopy(pmd_1["gen"]["1"])
    #
    #         gen["index"] = i+1
    #         gen["cost"] *= 0.01
    #         gen["configuration"] = load["configuration"]
    #         gen["pmax"] = (-load["pd"])/10
    #         gen["pmin"] *= 0
    #         gen["qmin"] = -abs.(gen["pmax"])/10
    #         gen["qmax"] = abs.(gen["pmax"])/10
    #         gen["gen_bus"] = load["load_bus"]
    #         gen["model"] = 2
    #     end
    #
    #     pm_ivr  = instantiate_mc_model(pmd, IVRUPowerModel, PMD.build_mc_opf)
    #     sol_ivr = optimize_model!(pm_ivr, optimizer=ipopt_solver)
    #     @assert(sol_1["termination_status"]==LOCALLY_SOLVED)
    #
    #     pm_acr  = instantiate_mc_model(pmd, ACRUPowerModel, PMD.build_mc_opf)
    #     sol_acr = optimize_model!(pm_acr, optimizer=ipopt_solver)
    #     @assert(sol_2["termination_status"]==LOCALLY_SOLVED)
    #
    # end
end
