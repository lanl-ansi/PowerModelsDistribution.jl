# tests of PowerModels generic functions

@testset "test PowerModels generic functions" begin
    @testset "build PMs.ref" begin
        data = TPPMs.parse_file("../test/data/matlab/case5_i_r_b.m")

        ref = PowerModels.build_ref(data)

        #TODO add some tests here
    end
end
