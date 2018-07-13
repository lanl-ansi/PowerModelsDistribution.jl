@testset "test idempotent units transformations" begin
    @testset "5-bus case" begin
        data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m")
        data_base = deepcopy(data)

        PMs.make_mixed_units(data)
        PMs.make_per_unit(data)
        @test InfrastructureModels.compare_dict(data, data_base)
    end
end


@testset "test angle wrapper functions" begin
    @testset "test wraptopi" begin
        wrappedradians = ThreePhasePowerModels.wraptopi([0, pi/2, pi, 3pi/2, 2pi])
        @test isapprox(wrappedradians, [0, pi/2, -pi, -pi/2, 0]; atol=1e-12)
    end

    @testset "test wrapto180" begin
        wrappeddegrees = ThreePhasePowerModels.wrapto180([0, 90, 180, 270, 360])
        @test isapprox(wrappeddegrees, [0, 90, -180, -90, 0]; atol=1e-12)
    end
end
