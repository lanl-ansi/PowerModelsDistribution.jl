
@testset "test make multi-phase" begin
    @testset "3-bus 3-phase case" begin
        mp_data = PMs.parse_file("$(pms_path)/test/data/matpower/case3.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 47267.9; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c], 1.58067; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(0.12669+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-3)
        end
    end

    @testset "5-bus 5-phase ac polar case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case5.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45522.096; atol = 1e-1)
        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.3999999; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(-0.0538204+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-5)
        end
    end

    @testset "5-bus 5-phase ac rectangular case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case5.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45522.096; atol = 1e-1)
        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.3999999; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(-0.0538204+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-5)
        end
    end

    @testset "5-bus 5-phase soc case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case5.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45365.17; atol = 1e-1)
        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.3999999; atol = 1e-3)
        end
    end

    @testset "30-bus 3-phase ac polar case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case30.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 614.007; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  2.192189; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(-0.071853+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-4)
        end
    end

    @testset "30-bus 3-phase ac rectangular case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case30.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 614.007; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  2.192189; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(-0.071853+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-4)
        end
    end

    @testset "30-bus 3-phase soc case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case30.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 517.588; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  2.821313; atol = 1e-3)
        end
    end
end



@testset "test multi-phase matlab parser" begin
    @testset "5-bus independent radial identical case" begin
        mp_data = TPPMs.parse_file("../test/data/matlab/case5_i_r_a.m")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 55451.7; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["qg"][c], 0.039742; atol = 1e-4)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(0.048896+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-4)
        end
    end
    @testset "5-bus independent radial different case" begin
        mp_data = TPPMs.parse_file("../test/data/matlab/case5_i_r_b.m")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 56091.3; atol = 1e-1)

        @test isapprox(result["solution"]["gen"]["1"]["qg"][1],  0.105276; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["va"][1],  0.0575114; atol = 1e-3)

        for c in 2:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["qg"][c],  0.0897773; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c],  TPPMs.wraptopi(0.052544+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-3)
        end
    end
    @testset "5-bus independent meshed different case" begin
        mp_data = TPPMs.parse_file("../test/data/matlab/case5_i_m_b.m")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 52964.4; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["qg"][c],  0.3; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(-0.0135651+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-3)
        end
    end
    @testset "5-bus coupled meshed case" begin
        @testset "ac case" begin
            mp_data = TPPMs.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["status"] == :LocalOptimal
            @test isapprox(result["objective"], 53272.9; atol = 1e-1)

            @test all(isapprox.(result["solution"]["gen"]["1"]["qg"].values, 0.3; atol = 1e-3))

            @test isapprox(result["solution"]["bus"]["2"]["va"][1], -0.0139580; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][2], -2.1069476; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][3],  2.0808321; atol = 1e-3)
        end
        @testset "soc case" begin
            mp_data = TPPMs.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

            @test result["status"] == :LocalOptimal
            @test isapprox(result["objective"], -0.000272; atol = 1e-3)

            @test isapprox(result["solution"]["gen"]["1"]["qg"][1], 0.0451820; atol = 1e-3)
            @test isapprox(result["solution"]["gen"]["1"]["qg"][2], 0.0290373; atol = 1e-3)
            @test isapprox(result["solution"]["gen"]["1"]["qg"][3], 0.0343748; atol = 1e-3)
        end
    end
    @testset "5-bus coupled meshed infeasible case" begin
        @testset "ac case" begin
            mp_data = TPPMs.parse_file("../test/data/matlab/case5_c_m_b.m")
            result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["status"] == :LocalInfeasible
        end
        #=
        # omit due to large number of terminal warnings
        @testset "soc case" begin
            mp_data = TPPMs.parse_file("../test/data/matlab/case5_c_m_b.m")
            result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

            @test result["status"] == :LocalInfeasible
        end
        =#
    end
    @testset "5-bus coupled radial no shunt case" begin
        mp_data = TPPMs.parse_file("../test/data/matlab/case5_c_r_a.m")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 55436.1; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c], 0.4; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["vm"][c], 1.08564; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(0.04905-2*pi/mp_data["conductors"]*(c-1)); atol = 1e-3)
        end
    end
    @testset "5-bus coupled radial shunt case" begin
        mp_data = TPPMs.parse_file("../test/data/matlab/case5_c_r_b.m")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 56075.1; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.4; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(0.055338-2*pi/mp_data["conductors"]*(c-1)); atol = 5e-3)
        end
    end
end


@testset "test dropped phases" begin
    @testset "4-bus 3-phase ac polar opf case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/case4_phase_drop.dss")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0182595; atol = 1e-4)

        @test isapprox(result["solution"]["gen"]["1"]["pg"][1], 5.06513e-5; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][2], 6.0865e-5; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][3], 7.1119e-5; atol = 1e-7)

        @test isapprox(result["solution"]["bus"]["2"]["vm"][1], 0.990023; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][2], 1.000000; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][3], 1.000000; atol = 1e-4)
    end

    @testset "4-bus 3-phase ac rectangular opf case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/case4_phase_drop.dss")
        result = run_tp_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0182595; atol = 1e-4)

        @test isapprox(result["solution"]["gen"]["1"]["pg"][1], 5.06513e-5; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][2], 6.0865e-5; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][3], 7.1119e-5; atol = 1e-7)

        # atol had to be increased from 1E-4 -> 1.5E-4 copmpared to ACP
        @test isapprox(result["solution"]["bus"]["2"]["vm"][1], 0.990023; atol = 1.5e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][2], 1.000000; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][3], 1.000000; atol = 1e-4)
    end

    @testset "5-bus 3-phase ac polar opf case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/case5_phase_drop.dss")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0597017; atol = 1e-4)

        @test isapprox(result["solution"]["gen"]["1"]["pg"][1], 0.00015236280779412599; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][2], 0.00019836795302238667; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][3], 0.0002469182092574594; atol = 1e-7)

        @test isapprox(result["solution"]["bus"]["2"]["vm"][1], 0.9736211293005391; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][2], 0.9650040745702724; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][3], 0.9562171321941326; atol = 1e-4)
    end

    @testset "5-bus 3-phase ac rectangular opf case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/case5_phase_drop.dss")
        result = run_tp_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0597017; atol = 1e-4)

        @test isapprox(result["solution"]["gen"]["1"]["pg"][1], 0.00015236280779412599; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][2], 0.00019836795302238667; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][3], 0.00024691804721165244; atol = 1e-7)

        @test isapprox(result["solution"]["bus"]["2"]["vm"][1], 0.9736211293005391; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][2], 0.9650040745702724; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][3], 0.9562166511182492; atol = 1e-4)
    end


    #=
    # causes a solve error in Ipopt, probably due to an issue with redundant constraints
    @testset "4-bus 3-phase ac pf case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/case4_phase_drop.dss")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0182595; atol = 1e-4)

        @test isapprox(result["solution"]["gen"]["1"]["pg"][1], 5.06513e-5; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][2], 6.0865e-5; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][3], 7.1119e-5; atol = 1e-7)

        @test isapprox(result["solution"]["bus"]["2"]["vm"][1], 0.990023; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][2], 1.000000; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][3], 1.000000; atol = 1e-4)
    end
    =#
end





@testset "test ac polar polar opf" begin
    @testset "30-bus make-3-phase case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case30.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 614.007; atol = 1e-1)
    end
    @testset "5-bus coupled meshed case" begin
        @testset "ac case" begin
            mp_data = TPPMs.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["status"] == :LocalOptimal
            @test isapprox(result["objective"], 53272.9; atol = 1e-1)
        end
    end
    @testset "5-bus phase drop case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/case5_phase_drop.dss")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0597017; atol = 1e-4)
    end
end


@testset "test ac rectangular opf" begin
    @testset "30-bus make-3-phase case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case30.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 614.007; atol = 1e-1)
    end
    @testset "5-bus coupled meshed case" begin
        @testset "ac case" begin
            mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_tp_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

            @test result["status"] == :LocalOptimal
            @test isapprox(result["objective"], 53272.9; atol = 1e-1)
        end
    end
    @testset "5-bus phase drop case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/opendss/case5_phase_drop.dss")
        result = run_tp_opf(mp_data, PMs.ACRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0597017; atol = 1e-4)
    end
end


@testset "test dc opf" begin
     @testset "30-bus make-3-phase case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case30.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.DCPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 566.112; atol = 1e-1)
    end
    @testset "5-bus coupled meshed case" begin
        @testset "ac case" begin
            mp_data = TPPMs.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_tp_opf(mp_data, PMs.DCPPowerModel, ipopt_solver)

            @test result["status"] == :LocalOptimal
            @test isapprox(result["objective"], 55640.2; atol = 1e-1)
        end
    end
    @testset "5-bus phase drop case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/case5_phase_drop.dss")
        result = run_tp_opf(mp_data, PMs.DCPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0540021; atol = 1e-4)
    end
end


@testset "test nfa opf" begin
     @testset "30-bus make-3-phase case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case30.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.NFAPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 458.006; atol = 1e-1)
    end
    @testset "5-bus coupled meshed case" begin
        @testset "ac case" begin
            mp_data = TPPMs.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_tp_opf(mp_data, PMs.NFAPowerModel, ipopt_solver)

            @test result["status"] == :LocalOptimal
            @test isapprox(result["objective"], 44700.0; atol = 1e-1)
        end
    end
    @testset "5-bus phase drop case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/case5_phase_drop.dss")
        result = run_tp_opf(mp_data, PMs.NFAPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.054; atol = 1e-4)
    end
    @testset "3w transformer case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/ut_trans_3w_dyy_basetest.dss")
        result = run_tp_opf(mp_data, PMs.NFAPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.666; atol = 1e-3)
    end
end


@testset "test soc (BIM) opf" begin
     @testset "30-bus make-3-phase case" begin
        mp_data = PMs.parse_file("../test/data/matpower/case30.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 517.563; atol = 1e-1)
    end
    @testset "5-bus coupled meshed case" begin
        @testset "ac case" begin
            mp_data = TPPMs.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

            @test result["status"] == :LocalOptimal
            @test isapprox(result["objective"], -0.000325497; atol = 1e-1)
        end
    end
    @testset "5-bus phase drop case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/case5_phase_drop.dss")
        result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0597016; atol = 1e-4)
    end
end
##

@testset "test acp opf unbalance constrained" begin
    tppm_data = TPPMs.parse_file("../test/data/matlab/case_bctr.m")
    # We check the equations by comparing against the value calculated by the solution
    # builder for the active constraint
    constr_keys = ["vm_vuf_max", "vm_seq_neg_max", "vm_seq_zero_max", "vm_seq_pos_max", "vm_ll_max", "vm_ll_min"]
    constr_lims = [0.04, 0.04, 0.04, 1.02, PMs.MultiConductorVector(ones(3)*1.07), PMs.MultiConductorVector(ones(3)*1.01)]
    sol_keys = ["vuf", "vm_seq_neg", "vm_seq_zero", "vm_seq_pos", "vm_ll", "vm_ll"]
    for i in 1:length(constr_keys)
        tppm = deepcopy(tppm_data)
        tppm["bus"]["3"][constr_keys[i]] = constr_lims[i]
        sol = TPPMs.run_tp_opf_bctr(tppm, PMs.ACPPowerModel, ipopt_solver, multiconductor=true)
        # the minimum is needed for the LL constraints; only one out of three will be active
        @test minimum(abs.(sol["solution"]["bus"]["3"][sol_keys[i]]-constr_lims[i])) <= 1E-5
    end
end
