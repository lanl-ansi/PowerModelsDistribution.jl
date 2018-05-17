@testset "test_multiphase ac ots" begin
    @testset "5-bus independent radial w/ shunts" begin
        result = run_tp_ots("../test/data/matlab/case5_i_r_b.m", PMs.ACPPowerModel, juniper_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 56091.3; atol=1e-1)
    end

    @testset "3-bus 3-phase case" begin
        mp_data = PMs.parse_file("../test/data/matlab/case3.m")
        PMs.make_multiphase(mp_data, 3)
        result = run_tp_ots(mp_data, PMs.ACPPowerModel, juniper_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 17720.6; atol=1e-1)
    end

    # @testset "5-bus 3-phase case" begin
    #     mp_data = PMs.parse_file("../test/data/matlab/case5.m")
    #     PMs.make_multiphase(mp_data, 3)
    #     result = run_tp_ots(mp_data, PMs.ACPPowerModel, juniper_solver)

    #     @test result["status"] == :LocalOptimal
    #     @test isapprox(result["objective"], 45522.1; atol=1e-1)
    # end

end

@testset "test multiphase dc ots" begin
    @testset "5-bus 3-phase case" begin
        @testset "5-bus independent radial w/o shunts" begin
            result = run_tp_ots("../test/data/matlab/case5_i_r_a.m", PMs.DCPPowerModel, pajarito_solver)

            @test result["status"] == :Optimal
            @test isapprox(result["objective"], 54870.0; atol=1e-1)
        end

        @testset "5-bus independent radial w/ shunts" begin
            result = run_tp_ots("../test/data/matlab/case5_i_r_b.m", PMs.DCPPowerModel, pajarito_solver)

            @test result["status"] == :Optimal
            @test isapprox(result["objective"], 55410.0; atol=1e-1)
        end

        @testset "3-bus 3-phase case" begin
            mp_data = PMs.parse_file("../test/data/matlab/case3.m")
            PMs.make_multiphase(mp_data, 3)
            result = run_tp_ots(mp_data, PMs.DCPPowerModel, pajarito_solver)

            @test result["status"] == :Optimal
            @test isapprox(result["objective"], 17346.1; atol=1e-1)
        end

        @testset "5-bus 3-phase case" begin
            mp_data = PMs.parse_file("../test/data/matlab/case5.m")
            PMs.make_multiphase(mp_data, 3)
            result = run_tp_ots(mp_data, PMs.DCPPowerModel, pajarito_solver)

            @test result["status"] == :Optimal
            @test isapprox(result["objective"], 44973.7; atol=1e-1)
        end

    end
end