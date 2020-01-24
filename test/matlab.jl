@info "running matlab parser tests"

@testset "test matlab data parser" begin
    @testset "5-bus matlab minimal data parse" begin
        data = PMD.parse_file("../test/data/matlab/case5_i_r_a.m")

        @test length(data["bus"]) == 5
        @test length(data["load"]) == 3
        @test length(data["shunt"]) == 0
        @test length(data["gen"]) == 5
        @test length(data["branch"]) == 4

        @test isa(data["bus"]["1"]["vm"], MultiConductorVector{Float64})
        @test isa(data["load"]["1"]["pd"], MultiConductorVector{Float64})
        @test isa(data["gen"]["1"]["pg"], MultiConductorVector{Float64})
        @test isa(data["branch"]["1"]["b_fr"], PMs.MultiConductorMatrix{Float64})
        @test isa(data["branch"]["1"]["br_x"], PMs.MultiConductorMatrix{Float64})

        @test haskey(data["bus"]["1"], "bus_name")
    end

    @testset "5-bus matlab shunt data parse" begin
        data = PMD.parse_file("../test/data/matlab/case5_i_r_b.m")

        @test length(data["bus"]) == 5
        @test length(data["load"]) == 3
        @test length(data["shunt"]) == 2
        @test length(data["gen"]) == 5
        @test length(data["branch"]) == 4

        @test isa(data["shunt"]["1"]["bs"], MultiConductorVector{Float64})
    end
end
