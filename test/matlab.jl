TESTLOG = getlogger(PowerModels)

@testset "test matlab data parser" begin
    @testset "5-bus minimal data" begin
        data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_a.m")

        @test length(data["bus"]) == 5
        @test length(data["load"]) == 3
        @test length(data["shunt"]) == 0
        @test length(data["gen"]) == 5
        @test length(data["branch"]) == 4

        @test isa(data["bus"]["1"]["vm"], PMs.MultiConductorVector{Float64})
        @test isa(data["load"]["1"]["pd"], PMs.MultiConductorVector{Float64})
        @test isa(data["gen"]["1"]["pg"], PMs.MultiConductorVector{Float64})
        @test isa(data["branch"]["1"]["b_fr"], PMs.MultiConductorVector{Float64})
        @test isa(data["branch"]["1"]["br_x"], PMs.MultiConductorMatrix{Float64})

        @test haskey(data["bus"]["1"], "bus_name")
    end

    @testset "5-bus shunt data" begin
        data = ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m")

        @test length(data["bus"]) == 5
        @test length(data["load"]) == 3
        @test length(data["shunt"]) == 2
        @test length(data["gen"]) == 5
        @test length(data["branch"]) == 4

        @test isa(data["shunt"]["1"]["bs"], PMs.MultiConductorVector{Float64})
    end

    #=
    @testset "version warning" begin
        setlevel!(TESTLOG, "warn")

        @test_warn(TESTLOG, "matlab source data has unrecognized version 0.0.0, cannot translate to version 1.0.0, parse may be invalid",
                   ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m"))

        @test_warn(TESTLOG, "No version number found, file may not be compatible with parser.",
                   ThreePhasePowerModels.parse_file("../test/data/matlab/case5_i_r_b.m"))

        setlevel!(TESTLOG, "error")
    end
    =#
end
