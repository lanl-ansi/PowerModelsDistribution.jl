
@testset "test make multi-phase" begin
    @testset "3-bus 3-phase case" begin
        mp_data = PMs.parse_file("../../PowerModels/test/data/matpower/case3.m")
        PMs.make_multiphase(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 47267.9; atol = 1e-1)

        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph], 1.58067; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], 0.12669+2*pi/mp_data["phases"]*(ph-1); atol = 1e-3)
        end
    end
    @testset "5-bus 5-phase ac case" begin
        mp_data = PMs.parse_file("../test/data/matlab/case5.m")
        PMs.make_multiphase(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45522.096; atol = 1e-1)
        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph],  0.3999999; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], -0.0538204+2*pi/mp_data["phases"]*(ph-1); atol = 1e-5)
        end
    end

    @testset "5-bus 5-phase soc case" begin
        mp_data = PMs.parse_file("../test/data/matlab/case5.m")
        PMs.make_multiphase(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45365.17; atol = 1e-1)
        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph],  0.3999999; atol = 1e-3)
        end
    end

    @testset "30-bus 3-phase ac case" begin
        mp_data = PMs.parse_file("../test/data/matlab/case30.m")
        PMs.make_multiphase(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 614.007; atol = 1e-1)

        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph],  2.192189; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], -0.071853+2*pi/mp_data["phases"]*(ph-1); atol = 1e-4)
        end
    end

    @testset "30-bus 3-phase soc case" begin
        mp_data = PMs.parse_file("../test/data/matlab/case30.m")
        PMs.make_multiphase(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 517.588; atol = 1e-1)

        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph],  2.821313; atol = 1e-3)
        end
    end
end



@testset "test multi-phase parser" begin
    @testset "5-bus independent radial identical case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_a.m")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 55451.7; atol = 1e-1)

        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["qg"][ph], 0.039742; atol = 1e-4)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], 0.048896+2*pi/mp_data["phases"]*(ph-1); atol = 1e-4)
        end
    end
    @testset "5-bus independent radial different case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 56091.3; atol = 1e-1)

        @test isapprox(result["solution"]["gen"]["1"]["qg"][1],  0.105276; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["va"][1],  0.0575114; atol = 1e-3)

        for ph in 2:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["qg"][ph],  0.0897773; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph],  0.052544+2*pi/mp_data["phases"]*(ph-1); atol = 1e-3)
        end
    end
    @testset "5-bus independent meshed different case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_m_b.m")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 52964.4; atol = 1e-1)

        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["qg"][ph],  0.3; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], -0.0135651+2*pi/mp_data["phases"]*(ph-1); atol = 1e-3)
        end
    end
    @testset "5-bus coupled meshed case" begin
        @testset "ac case" begin
            mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

            @test result["status"] == :LocalOptimal
            @test isapprox(result["objective"], 53350.44; atol = 1e-1)

            for ph in 1:mp_data["phases"]
                @test isapprox(result["solution"]["gen"]["1"]["qg"][ph],  0.29999999; atol = 1e-3)
            end
            @test isapprox(result["solution"]["bus"]["2"]["va"][1], -0.0135573; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][2], -0.0123172+2*pi/mp_data["phases"]; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][3], -0.0136547+4*pi/mp_data["phases"]; atol = 1e-3)
        end
        @testset "soc case" begin
            mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_m_a.m")
            result = run_tp_opf(mp_data, PMs.SOCWRPowerModel, ipopt_solver)

            @test result["status"] == :LocalOptimal
            @test isapprox(result["objective"], 21553.9; atol = 1e-1)

            @test isapprox(result["solution"]["gen"]["1"]["qg"][1], -0.037881; atol = 1e-3)
            @test isapprox(result["solution"]["gen"]["1"]["qg"][2], -0.299998; atol = 1e-3)
            @test isapprox(result["solution"]["gen"]["1"]["qg"][3], -0.075929; atol = 1e-3)
        end
    end
    @testset "5-bus coupled radial no shunt case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_r_a.m")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 55448.1; atol = 1e-1)

        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph], 0.4; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["vm"][ph], 1.08564; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], 0.04905+2*pi/mp_data["phases"]*(ph-1); atol = 1e-3)
        end
    end
    @testset "5-bus coupled radial shunt case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_c_r_b.m")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 56087.4; atol = 1e-1)

        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph],  0.4; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], 0.055338+2*pi/mp_data["phases"]*(ph-1); atol = 5e-3)
        end
    end
end