
@testset "test make multi-phase" begin
    @testset "3-bus 3-phase case" begin
        mp_data = PMs.parse_file("../test/data/matlab/case3.m")
        PMs.make_multiphase(mp_data, 3)
        ThreePhasePowerModels.shift_phase_angles!(mp_data)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 17720.6; atol = 1e-1)

        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph], 1.58067; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], 0.12669+2*pi/mp_data["phases"]*(ph-1); atol = 1e-3)
        end
    end
    @testset "5-bus 5-phase case" begin
        mp_data = PMs.parse_file("../test/data/matlab/case5.m")
        PMs.make_multiphase(mp_data, 3)
        ThreePhasePowerModels.shift_phase_angles!(mp_data)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45522.096; atol = 1e-1)
        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph],  0.3999999; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], -0.0538204+2*pi/mp_data["phases"]*(ph-1); atol = 1e-5)
        end
    end
    @testset "30-bus 3-phase case" begin
        mp_data = PMs.parse_file("../test/data/matlab/case30.m")
        PMs.make_multiphase(mp_data, 3)
        ThreePhasePowerModels.shift_phase_angles!(mp_data)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 614.007; atol = 1e-1)

        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph],  2.192189; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], -0.071853+2*pi/mp_data["phases"]*(ph-1); atol = 1e-4)
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
        @test isapprox(result["solution"]["gen"]["1"]["qg"][2],  0.0897773; atol = 1e-3)
        @test isapprox(result["solution"]["gen"]["1"]["qg"][3],  0.0897773; atol = 1e-3)

        @test isapprox(result["solution"]["bus"]["2"]["va"][1],  0.0575114; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["va"][2],  0.052544+2*pi/mp_data["phases"]; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["va"][3],  0.052544+4*pi/mp_data["phases"]; atol = 1e-3)
    end
    @testset "5-bus independent meshed different case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_m_b.m")
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 52964.4; atol = 1e-1)

        @test isapprox(result["solution"]["gen"]["1"]["qg"][1],  0.3; atol = 1e-3)
        @test isapprox(result["solution"]["gen"]["1"]["qg"][2],  0.3; atol = 1e-3)
        @test isapprox(result["solution"]["gen"]["1"]["qg"][3],  0.3; atol = 1e-3)

        @test isapprox(result["solution"]["bus"]["2"]["va"][1], -0.0135651; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["va"][2], -0.0135651+2*pi/mp_data["phases"]; atol = 1e-3)
        @test isapprox(result["solution"]["bus"]["2"]["va"][3], -0.0135651+4*pi/mp_data["phases"]; atol = 1e-3)
    end
end