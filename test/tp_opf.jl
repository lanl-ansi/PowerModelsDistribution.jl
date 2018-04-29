
@testset "test ac api" begin
    @testset "3-bus 3-phase case" begin
        mp_data = PMs.parse_file("../../PowerModels/test/data/matpower/case3.m")
        PMs.make_multiphase(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 17720.6; atol = 1e-1)

        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph], 1.58067; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], 0.12669; atol = 1e-3)
        end
    end
    @testset "5-bus 5-phase case" begin
        mp_data = PMs.parse_file("../../PowerModels/test/data/matpower/case5.m")
        PMs.make_multiphase(mp_data, 5)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 91345.5; atol = 1e-1)
        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph],  0.4; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], -0.00103692; atol = 1e-5)
        end
    end
    @testset "30-bus 3-phase case" begin
        mp_data = PMs.parse_file("../../PowerModels/test/data/matpower/case30.m")
        PMs.make_multiphase(mp_data, 3)
        result = run_tp_opf(mp_data, PMs.ACPPowerModel, ipopt_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 614.905; atol = 1e-1)

        for ph in 1:mp_data["phases"]
            @test isapprox(result["solution"]["gen"]["1"]["pg"][ph],  2.18839; atol = 1e-3)
            @test isapprox(result["solution"]["bus"]["2"]["va"][ph], -0.071759; atol = 1e-4)
        end
    end
end
