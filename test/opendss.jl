TESTLOG = getlogger(PowerModels)

@testset "opendss parser" begin
    @testset "simple generator branch load" begin
        setlevel!(TESTLOG, "info")

        @test_warn(TESTLOG, "Not all OpenDSS features are supported, currently only minimal support for \
                      lines, loads, generators, and capacitors as shunts. Transformers and reactors \
                      as transformer branches are included, but value translation is not fully supported.",
                   ThreePhasePowerModels.parse_file("../test/data/opendss/simple.dss"))

        Memento.Test.@test_log(TESTLOG, "info", "Calling parse_dss on ../test/data/opendss/simple.dss",
                               ThreePhasePowerModels.parse_file("../test/data/opendss/simple.dss"))

        Memento.Test.@test_log(TESTLOG, "info", "Done parsing ../test/data/opendss/simple.dss",
                               ThreePhasePowerModels.parse_file("../test/data/opendss/simple.dss"))

        setlevel!(TESTLOG, "error")
    end
end
