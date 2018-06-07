TESTLOG = getlogger(PowerModels)

@testset "opendss parser" begin
    @testset "simple generator branch load" begin
        setlevel!(TESTLOG, "info")

        @test_warn(TESTLOG, "Not all OpenDSS features are supported, currently only minimal support for \
                      lines, loads, generators, and capacitors as shunts. Transformers and reactors \
                      as transformer branches are included, but value translation is not fully supported.",
                   ThreePhasePowerModels.parse_file("../test/data/opendss/test_simple.dss"))

        Memento.Test.@test_log(TESTLOG, "info", "Calling parse_dss on ../test/data/opendss/test_simple.dss",
                               ThreePhasePowerModels.parse_file("../test/data/opendss/test_simple.dss"))

        Memento.Test.@test_log(TESTLOG, "info", "Done parsing ../test/data/opendss/test_simple.dss",
                               ThreePhasePowerModels.parse_file("../test/data/opendss/test_simple.dss"))

        setlevel!(TESTLOG, "error")

        dss = ThreePhasePowerModels.parse_dss("../test/data/opendss/test_simple.dss")
        tppm = ThreePhasePowerModels.parse_file("../test/data/opendss/test_simple.dss")

        for (key, len) in zip(["bus", "gen", "branch", "load", "dcline"], [2, 1, 1, 1, 0])
            @test haskey(tppm, key)
            @test length(tppm[key]) == len
        end

        sol = ThreePhasePowerModels.run_tp_opf(tppm, PowerModels.SOCWRPowerModel, ipopt_solver)

        @test sol["status"] == :LocalOptimal
    end
    @testset "parser cases" begin
        setlevel!(TESTLOG, "info")

        @test_warn(TESTLOG, "Line 27 in \"test2_linecodes.dss\" contains an unsupported symbol, skipping",
                   ThreePhasePowerModels.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Command \"solve\" on line 57 in \"test2_master.dss\" is not supported, skipping.",
                   ThreePhasePowerModels.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Command \"show\" on line 59 in \"test2_master.dss\" is not supported, skipping.",
                   ThreePhasePowerModels.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "cannot parse \"x=(1.05 0.75 0.001 5 * - - 125 15.0 / sqr *)\" as Float64 in reactor, leaving as String.",
                   ThreePhasePowerModels.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "cannot parse \"tap=(0.01000  15 * 1 +)\" as Float64 in transformer, leaving as String.",
                   ThreePhasePowerModels.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "transformers are not yet supported, treating like non-transformer lines",
                   ThreePhasePowerModels.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "reactors as constant impedance elements is not yet supported, treating like line",
                   ThreePhasePowerModels.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.Test.@test_log(TESTLOG, "info", "`dss_data` has been reset with the \"clear\" command.",
                               ThreePhasePowerModels.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.Test.@test_log(TESTLOG, "info", "Redirecting to file \"test2_linecodes.dss\"",
                               ThreePhasePowerModels.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.Test.@test_log(TESTLOG, "info", "Compiling file \"test2_loadshape.dss\"",
                               ThreePhasePowerModels.parse_file("../test/data/opendss/test2_master.dss"))

        setlevel!(TESTLOG, "error")

        dss = ThreePhasePowerModels.parse_dss("../test/data/opendss/test2_master.dss")
        tppm = ThreePhasePowerModels.parse_file("../test/data/opendss/test2_master.dss")

        @test tppm["name"] == "TEST2"
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
