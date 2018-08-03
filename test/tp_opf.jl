
@testset "test make multi-phase" begin
    @testset "3-bus 3-phase case" begin
        mp_data = PMs.parse_file("../../PowerModels/test/data/matpower/case3.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 47267.9; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c], 1.58067; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(0.12669+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-3)
        end
    end
    @testset "5-bus 5-phase ac case" begin
        mp_data = PMs.parse_file("../test/data/matlab/case5.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45522.096; atol = 1e-1)
        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.3999999; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(-0.0538204+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-5)
        end
    end

    @testset "5-bus 5-phase soc case" begin
        mp_data = PMs.parse_file("../test/data/matlab/case5.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45365.17; atol = 1e-1)
        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.3999999; atol = 1e-3)
        end
    end

    @testset "30-bus 3-phase ac case" begin
        mp_data = PMs.parse_file("../test/data/matlab/case30.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 614.007; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  2.192189; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(-0.071853+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-4)
        end
    end

    @testset "30-bus 3-phase soc case" begin
        mp_data = PMs.parse_file("../test/data/matlab/case30.m")
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
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_a.m")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 55451.7; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["qg"][c], 0.039742; atol = 1e-4)
            @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(0.048896+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-4)
        end
    end
    @testset "5-bus independent radial different case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m")
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
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_m_b.m")
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
            mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["status"] == :LocalOptimal
            @test isapprox(result["objective"], 53273.28; atol = 1e-1)

            @test all(isapprox.(result["solution"]["gen"]["1"]["qg"].values, 0.3; atol = 1e-3))

            @test isapprox(result["solution"]["bus"]["2"]["va"][1], TPPMs.wraptopi(-0.0135573); atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][2], TPPMs.wraptopi(-0.0123172-2*pi/mp_data["conductors"]); atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][3], TPPMs.wraptopi(-0.0136547-4*pi/mp_data["conductors"]); atol = 1e-3)
        end
        @testset "soc case" begin
            mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

            @test result["status"] == :LocalOptimal
            @test isapprox(result["objective"], -0.000272; atol = 1e-3)

            @test isapprox(result["solution"]["gen"]["1"]["qg"][1], 0.0472219; atol = 1e-3)
            @test isapprox(result["solution"]["gen"]["1"]["qg"][2], 0.0325493; atol = 1e-3)
            @test isapprox(result["solution"]["gen"]["1"]["qg"][3], 0.0357746; atol = 1e-3)
        end
    end
    @testset "5-bus coupled meshed infeasible case" begin
        @testset "ac case" begin
            mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_b.m")
            result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["status"] == :LocalInfeasible
        end
        #=
        # omit due to large number of terminal warnings
        @testset "soc case" begin
            mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_b.m")
            result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

            @test result["status"] == :LocalInfeasible
        end
        =#
    end
    @testset "5-bus coupled radial no shunt case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_r_a.m")
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
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_r_b.m")
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
    @testset "4-bus 3-phase ac opf case" begin
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

    @testset "5-bus 3-phase ac opf case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/case5_phase_drop.dss")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0185406; atol = 1e-4)

        @test isapprox(result["solution"]["gen"]["1"]["pg"][1], 5.132752553077026e-5; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][2], 6.177139067714014e-5; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][3], 7.230678136903121e-5; atol = 1e-7)

        @test isapprox(result["solution"]["bus"]["2"]["vm"][1], 0.9899205823173715; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][2], 0.9879617649716861; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][3], 0.9859780454153353; atol = 1e-4)
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

    @testset "5-bus 3-phase ac pf case" begin
        mp_data = TPPMs.parse_file("../test/data/opendss/case5_phase_drop.dss")
        result = run_tp_pf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 0.0; atol = 1e-4)

        @test isapprox(result["solution"]["gen"]["1"]["pg"][1], 5.13385961984349e-5; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][2], 6.17829534652931e-5; atol = 1e-7)
        @test isapprox(result["solution"]["gen"]["1"]["pg"][3], 7.23305690855803e-5; atol = 1e-7)

        @test isapprox(result["solution"]["bus"]["2"]["vm"][1], 0.9898208875771091; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][2], 0.9878608864128781; atol = 1e-4)
        @test isapprox(result["solution"]["bus"]["2"]["vm"][3], 0.9858740381554243; atol = 1e-4)
    end

end
