@testset "test linearised distflow opf_bf" begin
    @testset "3-bus case" begin
        mp_data = PowerModels.parse_file("../../PowerModels/test/data/matpower/case3.m")
        PowerModels.make_multiconductor(mp_data, 3)
        result = run_tp_opf_bf(mp_data, LPLinUBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45500.2 ; atol = 1e0)
        @test isapprox(result["solution"]["bus"]["3"]["vm"].values, 0.974398*[1,1,1]; atol = 1e-3)
    end
    @testset "5-bus case" begin
        mp_data = PowerModels.parse_file("../test/data/matlab/case5.m")
        PowerModels.make_multiconductor(mp_data, 3)
        result = run_tp_opf_bf(mp_data, LPLinUBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 44880; atol = 1e0)
        @test isapprox(result["solution"]["bus"]["3"]["vm"].values, 0.911466*[1,1,1]; atol = 1e-3)

    end
    @testset "5-bus independent radial identical case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_a.m")
        result = run_tp_opf_bf(mp_data, LPLinUBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 54870.0; atol = 1e-1)
        @test isapprox(result["solution"]["bus"]["3"]["vm"].values, 1.01026*[1,1,1]; atol = 1e-3)
    end
    @testset "5-bus independent radial different case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m")
        result = run_tp_opf_bf(mp_data, LPLinUBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 55307.7; atol = 1e-1)
        @test isapprox(result["solution"]["bus"]["3"]["vm"].values, [0.930014, 0.930014, 0.930014]; atol = 1e-3)

    end
end

@testset "test linearised distflow opf_bf in diagonal matrix form" begin
    @testset "3-bus case" begin
        mp_data = PowerModels.parse_file("../../PowerModels/test/data/matpower/case3.m")
        PowerModels.make_multiconductor(mp_data, 3)
        result = run_tp_opf_bf(mp_data, LPdiagUBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45500.2 ; atol = 1e0)
        # @test isapprox(result["solution"]["bus"]["3"]["vm"].values, 0.974398*[1,1,1]; atol = 1e-3)
    end
    @testset "5-bus case" begin
        mp_data = PowerModels.parse_file("../test/data/matlab/case5.m")
        PowerModels.make_multiconductor(mp_data, 3)
        result = run_tp_opf_bf(mp_data, LPdiagUBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 44880; atol = 1e0)
        # @test isapprox(result["solution"]["bus"]["3"]["vm"].values, 0.911466*[1,1,1]; atol = 1e-3)

    end
    @testset "5-bus independent radial identical case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_a.m")
        result = run_tp_opf_bf(mp_data, LPdiagUBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 54870.0; atol = 1e-1)
        @test isapprox(result["solution"]["bus"]["3"]["vm"].values, 1.01026*[1,1,1]; atol = 1e-3)

    end
    @testset "5-bus independent radial different case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m")
        result = run_tp_opf_bf(mp_data, LPdiagUBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 55307.7; atol = 1e-1)
        @test isapprox(result["solution"]["bus"]["3"]["vm"].values, [0.930014, 0.930014, 0.930014]; atol = 1e-3)

    end
end

@testset "test linearised distflow opf_bf in full matrix form" begin
    @testset "3-bus case" begin
        mp_data = PowerModels.parse_file("../../PowerModels/test/data/matpower/case3.m")
        PowerModels.make_multiconductor(mp_data, 3)
        result = run_tp_opf_bf(mp_data, LPfullUBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45500.2 ; atol = 1e0)
        # @test isapprox(result["solution"]["bus"]["3"]["vm"].values, 0.974398*[1,1,1]; atol = 1e-3)

    end
    @testset "5-bus case" begin
        mp_data = PowerModels.parse_file("../test/data/matlab/case5.m")
        PowerModels.make_multiconductor(mp_data, 3)
        result = run_tp_opf_bf(mp_data, LPfullUBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 44880; atol = 1e0)
        # @test isapprox(result["solution"]["bus"]["3"]["vm"].values, 0.911466*[1,1,1]; atol = 1e-3)

    end
    @testset "5-bus independent radial identical case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_a.m")
        result = run_tp_opf_bf(mp_data, LPfullUBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 54870.0; atol = 1e-1)
        # @test isapprox(result["solution"]["bus"]["3"]["vm"].values, 1.01026*[1,1,1]; atol = 1e-3)

    end
    @testset "5-bus independent radial different case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m")
        result = run_tp_opf_bf(mp_data, LPfullUBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 55307.7; atol = 1e-1)
        @test isapprox(result["solution"]["bus"]["3"]["vm"].values, [0.930014, 0.930014, 0.930014]; atol = 1e-3)

    end
end

@testset "test sdp distflow opf_bf" begin
    @testset "5-bus independent radial identical case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_a.m")
        result = run_tp_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

        @test result["status"] == :Optimal
        # @test isapprox(result["objective"], 55451.7; atol = 1e-1)
        @test isapprox(result["objective"], 55453.76; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            # @test isapprox(result["solution"]["gen"]["1"]["qg"][c], 0.039742; atol = 1e-4)
            # @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(0.048896+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-4)
        end
    end
    @testset "5-bus independent radial different case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m")
        result = run_tp_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

        @test result["status"] == :Optimal
        # @test isapprox(result["objective"], 56091.3; atol = 1e-1)
        @test isapprox(result["objective"], 56094.23; atol = 1e-1)

        @test isapprox(result["solution"]["gen"]["1"]["qg"][1],  0.0928659; atol = 1e-3)
        # @test isapprox(result["solution"]["gen"]["1"]["qg"][1],  0.105276; atol = 1e-3)
        # @test isapprox(result["solution"]["bus"]["2"]["va"][1],  0.0575114; atol = 1e-3)

        for c in 2:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["qg"][c],  0.138; atol = 1e-3)
            # @test isapprox(result["solution"]["bus"]["2"]["va"][c],  TPPMs.wraptopi(0.052544+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-3)
        end
    end
    @testset "5-bus independent meshed different case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_m_b.m")
        result = run_tp_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

        @test result["status"] == :Optimal
        @test isapprox(result["objective"], 45326.5; atol = 1e-1)

        # for c in 1:mp_data["conductors"]
        # @test isapprox(result["solution"]["gen"]["1"]["qg"][c],  0.3; atol = 1e-3)
        @test isapprox(result["solution"]["gen"]["1"]["qg"].values,  [0.08645, 0.1295, 0.1301]; atol = 1e-3)
            # @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(-0.0135651+2*pi/mp_data["conductors"]*(1-c)); atol = 1e-3)
        end
    end
    @testset "5-bus coupled meshed case" begin
        @testset "ac case" begin
            mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_tp_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

            @test result["status"] == :Optimal
            # @test isapprox(result["objective"], 53273.28; atol = 1e-1)
            @test isapprox(result["objective"], 45552.66; atol = 1e-1)

            # @test all(isapprox.(result["solution"]["gen"]["1"]["qg"].values, 0.3; atol = 1e-3))
            @test all(isapprox.(result["solution"]["gen"]["1"]["qg"].values, [0.17912, 0.15860, 0.223877]; atol = 1e-3))


            # @test isapprox(result["solution"]["bus"]["2"]["va"][1], TPPMs.wraptopi(-0.0135573); atol = 1e-3)
            # @test isapprox(result["solution"]["bus"]["2"]["va"][2], TPPMs.wraptopi(-0.0123172-2*pi/mp_data["conductors"]); atol = 1e-3)
            # @test isapprox(result["solution"]["bus"]["2"]["va"][3], TPPMs.wraptopi(-0.0136547-4*pi/mp_data["conductors"]); atol = 1e-3)
        end
    end
    @testset "5-bus coupled meshed infeasible case" begin
        @testset "ac case" begin
            mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_b.m")
            result = run_tp_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

            @test result["status"] == :UnknownError
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
        result = run_tp_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

        @test result["status"] == :Optimal
        # @test isapprox(result["objective"], 55436.1; atol = 1e-1)
        @test isapprox(result["objective"], 55436.6; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c], 0.4; atol = 1e-3)
            # @test isapprox(result["solution"]["bus"]["2"]["vm"][c], 1.08564; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["vm"][c], 1.08619; atol = 1e-3)
            # @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(0.04905-2*pi/mp_data["conductors"]*(c-1)); atol = 1e-3)
        end
    end
    @testset "5-bus coupled radial shunt case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_r_b.m")
        result = run_tp_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

        @test result["status"] == :Optimal
        # @test isapprox(result["objective"], 56075.1; atol = 1e-1)
        @test isapprox(result["objective"], 56074.76; atol = 1e-1)

        for c in 1:mp_data["conductors"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.4; atol = 1e-3)
            # @test isapprox(result["solution"]["bus"]["2"]["va"][c], TPPMs.wraptopi(0.055338-2*pi/mp_data["conductors"]*(c-1)); atol = 5e-3)
        end
    end
