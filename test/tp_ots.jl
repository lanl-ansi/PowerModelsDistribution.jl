function check_br_status(sol)
    for (i,branch) in sol["branch"]
        for c in 1:length(branch["br_status"])
            @test isapprox(branch["br_status"][c], 0.0, rtol=1e-6) || isapprox(branch["br_status"][c], 1.0, rtol=1e-6)
        end
    end
end


@testset "test_multiphase ac ots" begin
    @testset "3-bus 3-phase case" begin
        mp_data = PMs.parse_file("$(pms_path)/test/data/matpower/case3.m")
        PMs.make_multiconductor(mp_data, 3)
        result = run_tp_ots(mp_data, PMs.ACPPowerModel, juniper_solver)

        @test result["status"] == :LocalOptimal
        @test isapprox(result["objective"], 17720.6; atol=1e-1)

        check_br_status(result["solution"])
    end
end

@testset "test multiphase dc ots" begin
    @testset "5-bus 3-phase case" begin
        @testset "5-bus independent radial w/o shunts" begin
            result = run_tp_ots("../test/data/matlab/case5_i_r_a.m", PMs.DCPPowerModel, pavito_solver)

            @test result["status"] == :Optimal
            @test isapprox(result["objective"], 54870.0; atol=1e-1)

            check_br_status(result["solution"])
        end

        @testset "5-bus independent radial w/ shunts" begin
            result = run_tp_ots("../test/data/matlab/case5_i_r_b.m", PMs.DCPPowerModel, pavito_solver)

            @test result["status"] == :Optimal
            @test isapprox(result["objective"], 55410.0; atol=1e-1)

            check_br_status(result["solution"])
        end

        @testset "3-bus 3-phase case" begin
            mp_data = PMs.parse_file("../../PowerModels/test/data/matpower/case3.m")
            PMs.make_multiconductor(mp_data, 3)
            result = run_tp_ots(mp_data, PMs.DCPPowerModel, pavito_solver)

            @test result["status"] == :Optimal
            @test isapprox(result["objective"], 17346.1; atol=1e-1)

            check_br_status(result["solution"])
        end

        @testset "5-bus 3-phase case" begin
            mp_data = PMs.parse_file("../test/data/matpower/case5.m")
            PMs.make_multiconductor(mp_data, 3)
            result = run_tp_ots(mp_data, PMs.DCPPowerModel, pavito_solver)

            @test result["status"] == :Optimal
            @test isapprox(result["objective"], 44973.7; atol=1e-1)

            check_br_status(result["solution"])
        end

    end
end