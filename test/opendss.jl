TESTLOG = getlogger(PowerModels)
TPPMs = ThreePhasePowerModels

@testset "opendss parser" begin
    @testset "reverse polish notation" begin
        # Examples from OpenDSS manual
        @test isapprox(TPPMs.parse_rpn("2 pi * 60 * .001 *"), 2 * pi * 60 * .001; atol=1e-12)
        @test isapprox(TPPMs.parse_rpn("(14.4 13.8 / sqr 300 *"), (14.4 / 13.8)^2 * 300; atol=1e-12)
        @test isapprox(TPPMs.parse_rpn("24.9 3 sqrt /"), 24.9 / sqrt(3); atol=1e-12)

        @test TPPMs.isa_rpn("2 pi * 60 * .001 *")
        @test !TPPMs.isa_rpn("[ 2 10 ]")

        setlevel!(TESTLOG, "warn")

        @test_warn(TESTLOG, "parse_rpn does not support \"rollup\", \"rolldn\", or \"swap\", leaving as String",
                   TPPMs.parse_rpn("1 2 swap atan2"))

        @test_warn(TESTLOG, "\" 1 2 + - \" is not valid Reverse Polish Notation, leaving as String",
                   TPPMs.parse_rpn(" 1 2 + - "))

        @test_warn(TESTLOG, "\"1 2 3 +\" is not valid Reverse Polish Notation, leaving as String",
                   TPPMs.parse_rpn("1 2 3 +"))

        setlevel!(TESTLOG, "error")
    end
    @testset "simple generator branch load" begin
        setlevel!(TESTLOG, "info")

        @test_warn(TESTLOG, "Not all OpenDSS features are supported, currently only minimal support for \
                      lines, loads, generators, and capacitors as shunts. Transformers and reactors \
                      as transformer branches are included, but value translation is not fully supported.",
                      TPPMs.parse_file("../test/data/opendss/test_simple.dss"))

        Memento.Test.@test_log(TESTLOG, "info", "Calling parse_dss on ../test/data/opendss/test_simple.dss",
                               TPPMs.parse_file("../test/data/opendss/test_simple.dss"))

        Memento.Test.@test_log(TESTLOG, "info", "Done parsing ../test/data/opendss/test_simple.dss",
                               TPPMs.parse_file("../test/data/opendss/test_simple.dss"))

        setlevel!(TESTLOG, "error")

        dss = TPPMs.parse_dss("../test/data/opendss/test_simple.dss")
        tppm = TPPMs.parse_file("../test/data/opendss/test_simple.dss")

        for (key, len) in zip(["bus", "gen", "branch", "load", "dcline"], [2, 1, 1, 1, 0])
            @test haskey(tppm, key)
            @test length(tppm[key]) == len
        end

        sol = TPPMs.run_tp_opf(tppm, PowerModels.SOCWRPowerModel, ipopt_solver)

        @test sol["status"] == :LocalOptimal
    end
    @testset "parser cases" begin
        setlevel!(TESTLOG, "info")

        @test_warn(TESTLOG, "Line 27 in \"test2_linecodes.dss\" contains an unsupported symbol, skipping",
                   TPPMs.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Command \"solve\" on line 57 in \"test2_master.dss\" is not supported, skipping.",
                   TPPMs.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Command \"show\" on line 59 in \"test2_master.dss\" is not supported, skipping.",
                   TPPMs.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "transformers are not yet supported, treating like non-transformer lines",
                   TPPMs.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "reactors as constant impedance elements is not yet supported, treating like line",
                   TPPMs.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.Test.@test_log(TESTLOG, "info", "`dss_data` has been reset with the \"clear\" command.",
                               TPPMs.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.Test.@test_log(TESTLOG, "info", "Redirecting to file \"test2_linecodes.dss\"",
                               TPPMs.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.Test.@test_log(TESTLOG, "info", "Compiling file \"test2_loadshape.dss\"",
                               TPPMs.parse_file("../test/data/opendss/test2_master.dss"))

        setlevel!(TESTLOG, "error")

        dss = TPPMs.parse_dss("../test/data/opendss/test2_master.dss")
        tppm = TPPMs.parse_file("../test/data/opendss/test2_master.dss")

        @test tppm["name"] == "test2"
        @test length(tppm) == 16
        @test length(dss) == 12

        for (key, len) in zip(["bus", "load", "shunt", "branch", "gen", "dcline"], [11, 3, 2, 10, 2, 0])
            @test haskey(tppm, key)
            @test length(tppm[key]) == len
        end

        for key in ["loadshape", "linecode", "buscoords", "options", "filename"]
            @test haskey(dss, key)
        end
    end
end
