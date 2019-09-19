@info "running opendss parser tests"

@testset "test opendss parser" begin
    @testset "reverse polish notation" begin
        @test isapprox(PMD._parse_rpn("2 pi * 60 * .001 *"), 2 * pi * 60 * .001; atol=1e-12)
        @test isapprox(PMD._parse_rpn("(14.4 13.8 / sqr 300 *"), (14.4 / 13.8)^2 * 300; atol=1e-12)
        @test isapprox(PMD._parse_rpn("24.9 3 sqrt /"), 24.9 / sqrt(3); atol=1e-12)
        @test all(isapprox.(PMD._parse_array(Float64, "(\"24.9 3 sqrt /\" \"10 2 *\")"), [24.9 / sqrt(3), 10 * 2.0]; atol=1e-12))

        @test PMD._isa_rpn("2 pi * 60 * .001 *")
        @test !PMD._isa_rpn("[ 2 10 ]")

        Memento.setlevel!(TESTLOG, "warn")

        @test_warn(TESTLOG, "_parse_rpn does not support \"rollup\", \"rolldn\", or \"swap\", leaving as String",
                   PMD._parse_rpn("1 2 swap atan2"))

        @test_warn(TESTLOG, "\" 1 2 + - \" is not valid Reverse Polish Notation, leaving as String",
                   PMD._parse_rpn(" 1 2 + - "))

        @test_warn(TESTLOG, "\"1 2 3 +\" is not valid Reverse Polish Notation, leaving as String",
                   PMD._parse_rpn("1 2 3 +"))

        Memento.setlevel!(TESTLOG, "error")
    end

    @testset "opendss parse load model errors" begin
        for load in 1:5
           dss = PMD.parse_dss("../test/data/opendss/loadparser_error.dss")
           dss["load"] = [dss["load"][load]]
           @test_throws(TESTLOG, ErrorException, PMD.parse_opendss(dss))
        end
    end

    @testset "opendss parse load model warnings" begin
        for model in [3, 4, 7, 8]
           dss = PMD.parse_dss("../test/data/opendss/loadparser_warn_model.dss")
           dss["load"] = [l for l in dss["load"] if l["name"]=="d1phm$model"]
           Memento.setlevel!(TESTLOG, "info")
           @test_warn(TESTLOG, ": load model $model not supported. Treating as model 1.", PMD.parse_opendss(dss))
           Memento.setlevel!(TESTLOG, "error")
        end
    end


    @testset "opendss parse generic warnings and errors" begin
        Memento.setlevel!(TESTLOG, "info")

        @test_warn(TESTLOG, "Not all OpenDSS features are supported, currently only minimal support for lines, loads, generators, and capacitors as shunts. Transformers and reactors as transformer branches are included, but value translation is not fully supported.",
                PMD.parse_file("../test/data/opendss/test_simple.dss"))

        Memento.TestUtils.@test_log(TESTLOG, "info", "Calling parse_dss on ../test/data/opendss/test_simple.dss",
                        PMD.parse_file("../test/data/opendss/test_simple.dss"))

        Memento.TestUtils.@test_log(TESTLOG, "info", "Done parsing ../test/data/opendss/test_simple.dss",
                        PMD.parse_file("../test/data/opendss/test_simple.dss"))

        @test_throws(TESTLOG, ErrorException,
                   PMD.parse_file("../test/data/opendss/test_simple2.dss"))

        @test_throws(TESTLOG, ErrorException,
                   PMD.parse_file("../test/data/opendss/test_simple2.dss"))

        @test_warn(TESTLOG, "Command \"solve\" on line 69 in \"test2_master.dss\" is not supported, skipping.",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Command \"show\" on line 71 in \"test2_master.dss\" is not supported, skipping.",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "reactors as constant impedance elements is not yet supported, treating like line",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "line.l1: like=something cannot be found",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Rg,Xg are not fully supported",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Could not find line \"something\"",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "The neutral impedance, (rneut and xneut properties), is ignored; the neutral (for wye and zig-zag windings) is connected directly to the ground.",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Only three-phase transformers are supported. The bus specification b7.1 is treated as b7 instead.",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Only three-phase transformers are supported. The bus specification b7.1 is treated as b7 instead.",
                  PMD.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.TestUtils.@test_log(TESTLOG, "info", "`dss_data` has been reset with the \"clear\" command.",
                               PMD.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.TestUtils.@test_log(TESTLOG, "info", "Redirecting to file \"test2_Linecodes.dss\"",
                               PMD.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.TestUtils.@test_log(TESTLOG, "info", "Compiling file \"test2_Loadshape.dss\"",
                               PMD.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.setlevel!(TESTLOG, "error")
    end

    pmd = PMD.parse_file("../test/data/opendss/test2_master.dss")

    len = 0.013516796
    rmatrix=PMD._parse_matrix(Float64, "[1.5000  |0.200000  1.50000  |0.250000  0.25000  2.00000  ]") * 3
    xmatrix=PMD._parse_matrix(Float64, "[1.0000  |0.500000  0.50000  |0.500000  0.50000  1.000000  ]") * 3
    cmatrix = PMD._parse_matrix(Float64, "[8.0000  |-2.00000  9.000000  |-1.75000  -2.50000  8.00000  ]") / 3

    @testset "opendss parse generic parser verification" begin
        dss = PMD.parse_dss("../test/data/opendss/test2_master.dss")

        @test pmd["name"] == "test2"

        @test length(pmd) == 21
        @test length(dss) == 12

        for (key, len) in zip(["bus", "load", "shunt", "branch", "gen", "dcline", "transformer"], [33, 4, 5, 27, 4, 0, 10])
            @test haskey(pmd, key)
            @test length(pmd[key]) == len
        end

        @test all(haskey(dss, key) for key in ["loadshape", "linecode", "buscoords", "options", "filename"])
    end

    @testset "opendss parse generic branch values verification" begin
        basekv_br3 = pmd["bus"][string(pmd["branch"]["3"]["f_bus"])]["base_kv"]
        @test all(isapprox.(pmd["branch"]["3"]["br_r"].values, rmatrix * len / basekv_br3^2 * pmd["baseMVA"]; atol=1e-6))
        @test all(isapprox.(pmd["branch"]["3"]["br_x"].values, xmatrix * len / basekv_br3^2 * pmd["baseMVA"]; atol=1e-6))
        @test all(isapprox.(pmd["branch"]["3"]["b_fr"].values, diag(basekv_br3^2 / pmd["baseMVA"] * 2.0 * pi * 60.0 * cmatrix * len / 1e9) / 2.0; atol=1e-6))

        for i in 6:7
            basekv_bri = pmd["bus"][string(pmd["branch"]["$i"]["f_bus"])]["base_kv"]
            @test all(isapprox.(diag(pmd["branch"]["$i"]["b_fr"].values), (3.4 * 2.0 + 1.6) / 3.0 * (basekv_bri^2 / pmd["baseMVA"] * 2.0 * pi * 60.0 / 1e9) / 2.0 / 3; atol=1e-6))
        end

        @test all(isapprox.(pmd["branch"]["1"]["br_r"].values.*(115/69)^2, diagm(0 => fill(6.3012e-8, 3)); atol=1e-12))
        @test all(isapprox.(pmd["branch"]["1"]["br_x"].values.*(115/69)^2, diagm(0 => fill(6.3012e-7, 3)); atol=1e-12))

        for k in ["qd", "pd"]
            @test all(isapprox.(pmd["load"]["4"][k].values, pmd["load"]["2"][k].values; atol=1e-12))
        end

        for k in ["gs", "bs"]
            @test all(isapprox.(pmd["shunt"]["2"][k].values, pmd["shunt"]["3"][k].values; atol=1e-12))
            @test all(isapprox.(pmd["shunt"]["4"][k].values, pmd["shunt"]["5"][k].values; atol=1e-12))
        end

        for k in keys(pmd["gen"]["3"])
            if !(k in ["gen_bus", "index", "name", "source_id", "active_phases"])
                if isa(pmd["gen"]["3"][k], PMs.MultiConductorValue)
                    @test all(isapprox.(pmd["gen"]["4"][k].values, pmd["gen"]["3"][k].values; atol=1e-12))
                else
                    @test all(isapprox.(pmd["gen"]["4"][k], pmd["gen"]["3"][k]; atol=1e-12))
                end
            end
        end

        for k in keys(pmd["branch"]["11"])
            if !(k in ["f_bus", "t_bus", "index", "name", "linecode", "source_id", "active_phases"])
                mult = 1.0
                if k in ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to"]
                    # compensation for the different voltage base
                    basekv_br3 = pmd["bus"][string(pmd["branch"]["3"]["f_bus"])]["base_kv"]
                    basekv_br8 = pmd["bus"][string(pmd["branch"]["8"]["f_bus"])]["base_kv"]
                    zmult = (basekv_br3/basekv_br8)^2
                    mult = (k in ["br_r", "br_x"]) ? zmult : 1/zmult
                end
                if isa(pmd["branch"]["11"][k], PMs.MultiConductorValue)
                    @test all(isapprox.(pmd["branch"]["3"][k].values.*mult, pmd["branch"]["8"][k].values; atol=1e-12))
                else
                    @test all(isapprox.(pmd["branch"]["3"][k].*mult, pmd["branch"]["8"][k]; atol=1e-12))
                end
            end
        end
    end

    @testset "opendss parse length units" begin
        @test pmd["branch"]["9"]["length"] == 1000.0 * len
        basekv_br9 = pmd["bus"][string(pmd["branch"]["9"]["f_bus"])]["base_kv"]
        @test all(isapprox.(pmd["branch"]["9"]["br_r"].values, rmatrix * len / basekv_br9^2 * pmd["baseMVA"]; atol=1e-6))
        @test all(isapprox.(pmd["branch"]["9"]["br_x"].values, xmatrix * len / basekv_br9^2 * pmd["baseMVA"]; atol=1e-6))
        @test all(isapprox.(pmd["branch"]["9"]["b_fr"].values, diag(basekv_br9^2 / pmd["baseMVA"] * 2.0 * pi * 60.0 * cmatrix * len / 1e9) / 2.0 / 3; atol=1e-6))
    end

    @testset "opendss parse switch length verify" begin
        @testset "branches with switches" begin
            @test pmd["branch"]["5"]["switch"]
            @test pmd["branch"]["5"]["length"] == 0.001
            @test all([pmd["branch"]["$i"]["switch"] == false for i in 1:4])
        end
    end

    @testset "opendss parse transformer parsing verify" begin
        dss_data = PMD.parse_dss("../test/data/opendss/test_transformer_formatting.dss")
        transformer = dss_data["transformer"][1]
        @test transformer["phases"] == "3"
        @test transformer["tap"] == "(0.00625 12 * 1 +)"
        @test transformer["tap_2"] == "1.5"
        @test transformer["%loadloss"] == "0.01"
        @test transformer["xhl"] == "0.02"
        @test transformer["kv_2"] == "12.47"
        @test transformer["conn_2"] == "wye"
        @test transformer["tap_3"] == "0.9"
        @test transformer["wdg_3"] == "3"

        PMD._apply_like!(dss_data["transformer"][3], dss_data, "transformer")
        @test dss_data["transformer"][3]["%loadloss"] == dss_data["transformer"][2]["%loadloss"]

        pmd_data = PMD.parse_file("../test/data/opendss/test_transformer_formatting.dss")
        @test all(all(pmd_data["transformer"]["$n"]["tm"] .== tm) for (n, tm) in zip(1:3, [1.075, 1.5, 0.9]))
    end

    @testset "opendss parse storage" begin
        pmd_storage = PMD.parse_file("../test/data/opendss/case3_balanced_battery.dss")
        for bat in values(pmd_storage["storage"])
            for key in ["energy", "storage_bus", "energy_rating", "charge_rating", "discharge_rating",
                        "charge_efficiency", "discharge_efficiency", "thermal_rating", "qmin", "qmax",
                        "r", "x", "p_loss", "q_loss", "status", "source_id", "active_phases"]
                @test haskey(bat, key)
                if key in ["x", "r", "qmin", "qmax", "thermal_rating"]
                    @test isa(bat[key], PowerModels.MultiConductorVector)
                end
            end
        end

        @test pmd_storage["storage"]["1"]["source_id"] == "storage.s1" && length(pmd_storage["storage"]["1"]["active_phases"]) == 3
    end

    @testset "opendss parse pvsystem" begin
        Memento.setlevel!(TESTLOG, "warn")
        @test_warn(TESTLOG, "Converting PVSystem \"pv1\" into generator with limits determined by OpenDSS property 'kVA'",
                PMD.parse_file("../test/data/opendss/case3_balanced_pv.dss"))
        Memento.setlevel!(TESTLOG, "error")
    end

    @testset "opendss parse verify source_id" begin
        @test pmd["shunt"]["1"]["source_id"] == "capacitor.c1" && length(pmd["shunt"]["1"]["active_phases"]) == 3
        @test pmd["shunt"]["4"]["source_id"] == "reactor.reactor3" && length(pmd["shunt"]["4"]["active_phases"]) == 3

        @test pmd["branch"]["1"]["source_id"] == "line.l1" && length(pmd["branch"]["1"]["active_phases"]) == 3
        @test pmd["transformer"]["1"]["source_id"] == "transformer.t4_1"  # winding indicated by _1
        @test pmd["branch"]["10"]["source_id"] == "reactor.reactor1" && length(pmd["branch"]["10"]["active_phases"]) == 3

        @test pmd["gen"]["1"]["source_id"] == "vsource.sourcegen" && length(pmd["gen"]["1"]["active_phases"]) == 3
        @test pmd["gen"]["2"]["source_id"] == "generator.g1" && length(pmd["gen"]["2"]["active_phases"]) == 3

        source_id = PMD._parse_dss_source_id(pmd["load"]["1"])
        @test source_id.dss_type == "load"
        @test source_id.dss_name == "ld1"
        @test all([n in source_id.active_phases for n in 1:2])

        @test all(haskey(component, "source_id") && haskey(component, "active_phases") for component_type in ["load", "branch", "shunt", "gen", "storage", "pvsystem"] for component in values(get(pmd, component_type, Dict())))
    end

    @testset "opendss parse verify order of properties on line" begin
        pmd1 = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
        pmd2 = PMD.parse_file("../test/data/opendss/case3_balanced_prop-order.dss")

        @test pmd1 == pmd2

        dss1 = PMD.parse_dss("../test/data/opendss/case3_balanced.dss")
        dss2 = PMD.parse_dss("../test/data/opendss/case3_balanced_prop-order.dss")

        @test dss1 != dss2
        @test all(a == b for (a, b) in zip(dss2["line"][1]["prop_order"],["name", "bus1", "bus2", "linecode", "rmatrix", "length"]))
        @test all(a == b for (a, b) in zip(dss2["line"][2]["prop_order"],["name", "bus1", "bus2", "like", "linecode", "length"]))
    end
end

@testset "test json parser" begin
    pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")

    io = PipeBuffer()
    JSON.print(io, pmd)
    pmd_json_file = PMD.parse_file(io)

    @test pmd == pmd_json_file
end
