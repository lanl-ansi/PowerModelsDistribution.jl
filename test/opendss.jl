TESTLOG = getlogger(PowerModels)
TPPMs = ThreePhasePowerModels

@testset "opendss parser" begin
    @testset "reverse polish notation" begin
        # Examples from OpenDSS manual
        @test isapprox(TPPMs.parse_rpn("2 pi * 60 * .001 *"), 2 * pi * 60 * .001; atol=1e-12)
        @test isapprox(TPPMs.parse_rpn("(14.4 13.8 / sqr 300 *"), (14.4 / 13.8)^2 * 300; atol=1e-12)
        @test isapprox(TPPMs.parse_rpn("24.9 3 sqrt /"), 24.9 / sqrt(3); atol=1e-12)
        @test all(isapprox.(TPPMs.parse_array(Float64, "(\"24.9 3 sqrt /\" \"10 2 *\")"), [24.9 / sqrt(3), 10 * 2.0]; atol=1e-12))

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

        for (key, len) in zip(["bus", "gen", "branch", "load", "dcline"], [2, 2, 1, 1, 0])
            @test haskey(tppm, key)
            @test length(tppm[key]) == len
        end

        sol = TPPMs.run_tp_opf(tppm, PowerModels.SOCWRPowerModel, ipopt_solver)

        @test sol["status"] == :LocalOptimal
    end

    @testset "parser cases" begin
        setlevel!(TESTLOG, "info")

        @test_warn(TESTLOG, "Command \"solve\" on line 59 in \"test2_master.dss\" is not supported, skipping.",
                   TPPMs.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Command \"show\" on line 61 in \"test2_master.dss\" is not supported, skipping.",
                   TPPMs.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "transformers are not yet supported, treating like non-transformer lines",
                   TPPMs.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "reactors as constant impedance elements is not yet supported, treating like line",
                   TPPMs.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Rg,Xg are not fully supported",
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

        for (key, len) in zip(["bus", "load", "shunt", "branch", "gen", "dcline"], [11, 3, 2, 12, 3, 0])
            @test haskey(tppm, key)
            @test length(tppm[key]) == len
        end

        for key in ["loadshape", "linecode", "buscoords", "options", "filename"]
            @test haskey(dss, key)
        end

        len = 0.013516796
        rmatrix=TPPMs.parse_matrix(Float64, "[1.5000  |0.200000  1.50000  |0.250000  0.25000  2.00000  ]")
        xmatrix=TPPMs.parse_matrix(Float64, "[1.0000  |0.500000  0.50000  |0.500000  0.50000  1.000000  ]")
        cmatrix = TPPMs.parse_matrix(Float64, "[8.0000  |-2.00000  9.000000  |-1.75000  -2.50000  8.00000  ]")

        @test all(isapprox.(tppm["branch"]["3"]["br_r"].values, rmatrix * len / tppm["basekv"]^2; atol=1e-6))
        @test all(isapprox.(tppm["branch"]["3"]["br_x"].values, xmatrix * len / tppm["basekv"]^2; atol=1e-6))
        @test all(isapprox.(tppm["branch"]["3"]["b_fr"].values, diag(tppm["basekv"]^2 / tppm["baseMVA"]^2 * 2.0 * pi * 60.0 * cmatrix * len / 1e9) / 2.0; atol=1e-6))

        for i in 6:7
            @test all(isapprox.(tppm["branch"]["$i"]["b_fr"].values, (3.4 * 2.0 + 1.6) / 3.0 * (tppm["basekv"]^2 / tppm["baseMVA"]^2 * 2.0 * pi * 60.0 / 1e9) / 2.0; atol=1e-6))
        end

        @test all(isapprox.(tppm["branch"]["1"]["br_r"].values, diagm(fill(2.1004e-10, 3)); atol=1e-12))
        @test all(isapprox.(tppm["branch"]["1"]["br_x"].values, diagm(fill(2.1004e-9, 3)); atol=1e-12))
    end

    @testset "3-bus balanced" begin
        tppm = TPPMs.parse_file("../test/data/opendss/case3_balanced.dss")
        @testset "OPF" begin
            @testset "ACP" begin
                sol = TPPMs.run_tp_opf(tppm, PMs.ACPPowerModel, ipopt_solver)

                @test sol["status"] == :LocalOptimal

                for bus in values(sol["solution"]["bus"])
                    @test all(isapprox.(bus["va"].values, [2 * pi / tppm["phases"] * (ph - 1) for ph in 1:tppm["phases"]]; atol=2e-2))
                    @test all(isapprox.(bus["vm"].values, 0.9959; atol=5e-1))
                end

                @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183456; atol=1e-3)
                @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923328; atol=1e-3)
            end
            @testset "SOC" begin
                sol = TPPMs.run_tp_opf(tppm, PMs.SOCWRPowerModel, ipopt_solver)

                @test sol["status"] == :LocalOptimal
            end
            @testset "LDF" begin
                sol = TPPMs.run_tp_opf_bf(tppm, PMs.SOCBFPowerModel, ipopt_solver)

                @test sol["status"] == :LocalOptimal
                @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183456; atol=1e-3)
                @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923328; atol=1e-3)
            end
        end
    end

    @testset "3-bus unbalanced" begin
    tppm = TPPMs.parse_file("../test/data/opendss/case3_unbalanced.dss")

        @testset "OPF" begin
            @testset "ACP" begin
                sol = TPPMs.run_tp_opf(tppm, PMs.ACPPowerModel, ipopt_solver)

                @test sol["status"] == :LocalOptimal

                for bus in values(sol["solution"]["bus"])
                    @test all(isapprox.(bus["va"].values, [2 * pi / tppm["phases"] * (ph - 1) for ph in 1:tppm["phases"]]; atol=2e-2))
                    @test all(isapprox.(bus["vm"].values, 0.9959; atol=5e-1))
                end

                @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214835; atol=1e-3)
                @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00932693; atol=1e-3)
            end
            @testset "SOC" begin
                sol = TPPMs.run_tp_opf(tppm, PMs.SOCWRPowerModel, ipopt_solver)

                @test sol["status"] == :LocalOptimal
            end
            @testset "LDF" begin
                sol = TPPMs.run_tp_opf_bf(tppm, PMs.SOCBFPowerModel, ipopt_solver)

                @test sol["status"] == :LocalOptimal
                @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214835; atol=1e-3)
                @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00932693; atol=1e-3)
            end
        end
    end
end
