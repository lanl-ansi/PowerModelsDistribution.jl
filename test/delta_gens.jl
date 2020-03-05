@info "running generator configuration tests"

@testset "test generator configuration" begin
    # This test checks the generators are connected properly by comparing them
    # to equivalent constant-power loads. This is achieved by fixing their bounds.
    @testset "ACP/ACR tests" begin
        pmd_1 = PMD.parse_file("$pmd_path/test/data/opendss/case3_delta_gens.dss")

        # convert to constant power loads
        for (_, load) in pmd_1["load"]
            load["model"] = "constant_power"
        end

        # create data model with equivalent generators
        pmd_2 = deepcopy(pmd_1)
        pmd_2["load"] = Dict()
        gen2load = Dict()
        for (i,(id, load)) in enumerate(pmd_1["load"])
            load = pmd_1["load"][id]
            gen = deepcopy(pmd_1["gen"]["1"])

            gen["index"] = i+1
            gen["cost"] *= 0
            gen["configuration"] = load["configuration"]
            gen["pmax"] = gen["pmin"] = -load["pd"]
            gen["qmin"] = gen["qmax"] = -load["qd"]
            gen["gen_bus"] = load["load_bus"]
            gen["model"] = 2
            gen2load["$(i+1)"] = id
            pmd_2["gen"]["$(i+1)"] = gen
        end

        # check ACP and ACR
        for (form, build_method) in zip(
                [PMs.ACPPowerModel, PMs.ACRPowerModel, PMs.IVRPowerModel],
                [PMD.build_mc_opf, PMD.build_mc_opf, PMD.build_mc_opf_iv]
            )
            pm_1  = PMs.instantiate_model(pmd_1, form, build_method, ref_extensions=[PMD.ref_add_arcs_trans!], multiconductor=true)
            sol_1 = PMs.optimize_model!(pm_1, optimizer=ipopt_solver)
            @assert(sol_1["termination_status"]==LOCALLY_SOLVED)

            pm_2  = PMs.instantiate_model(pmd_2, form, build_method, ref_extensions=[PMD.ref_add_arcs_trans!], multiconductor=true)
            sol_2 = PMs.optimize_model!(pm_2, optimizer=ipopt_solver)
            @assert(sol_2["termination_status"]==LOCALLY_SOLVED)

            # check that gens are equivalent to the loads
            for (gen, load) in gen2load
                pd_bus = sol_1["solution"]["load"][load]["pd"]
                qd_bus = sol_1["solution"]["load"][load]["qd"]
                pg_bus = sol_2["solution"]["gen"][gen]["pg"]
                qg_bus = sol_2["solution"]["gen"][gen]["qg"]
                @test(isapprox(pd_bus, -pg_bus, atol=1E-5))
                @test(isapprox(qd_bus, -qg_bus, atol=1E-5))
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
    #     pm_ivr  = PMs.instantiate_model(pmd, PMs.IVRPowerModel, PMD.build_mc_opf_iv, ref_extensions=[PMD.ref_add_arcs_trans!], multiconductor=true)
    #     sol_ivr = PMs.optimize_model!(pm_ivr, optimizer=ipopt_solver)
    #     @assert(sol_1["termination_status"]==LOCALLY_SOLVED)
    #
    #     pm_acr  = PMs.instantiate_model(pmd, PMs.ACRPowerModel, PMD.build_mc_opf, ref_extensions=[PMD.ref_add_arcs_trans!], multiconductor=true)
    #     sol_acr = PMs.optimize_model!(pm_acr, optimizer=ipopt_solver)
    #     @assert(sol_2["termination_status"]==LOCALLY_SOLVED)
    #
    # end
end
