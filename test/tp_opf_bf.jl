@testset "test soc distflow opf_bf" begin
    @testset "3-bus case" begin
        mp_data = PowerModels.parse_file("../../PowerModels/test/data/matpower/case3.m")
        PowerModels.make_multiconductor(mp_data, 3)
        result = run_tp_opf_bf(mp_data, PMs.SOCBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45895.9 ; atol = 1e0)
    end
    @testset "5-bus case" begin
        mp_data = PowerModels.parse_file("../test/data/matlab/case5.m")
        PowerModels.make_multiconductor(mp_data, 3)
        result = run_tp_opf_bf(mp_data, PMs.SOCBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 45069.6; atol = 1e0)
    end
    @testset "5-bus independent radial identical case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_a.m")
        result = run_tp_opf_bf(mp_data, PMs.SOCBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 55093.06; atol = 1e-1)
    end
    @testset "5-bus independent radial different case" begin
        mp_data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m")
        result = run_tp_opf_bf(mp_data, PMs.SOCBFPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 55636.2; atol = 1e-1)
    end
end
