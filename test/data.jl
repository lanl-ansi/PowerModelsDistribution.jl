@testset "test idempotent units transformations" begin
    @testset "5-bus case" begin
        data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m")
        data_base = deepcopy(data)

        PMs.make_mixed_units(data)
        PMs.make_per_unit(data)
        @test InfrastructureModels.compare_dict(data, data_base)
    end
end