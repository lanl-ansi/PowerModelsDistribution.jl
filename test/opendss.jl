@info "running opendss parser tests"

@testset "test opendss parser" begin

    @testset "loadshape parsing" begin
        dss = PMD.parse_dss("../test/data/opendss/loadshapes.dss")

        loadshapes = Dict{String,Any}()
        for (name, ls) in dss["loadshape"]
            loadshapes[name] = PMD._create_loadshape(name; PMD._to_kwargs(ls)...)
        end

        @test isapprox(loadshapes["1"]["interval"], 1.0/60)
        @test all(length(ls["pmult"]) == 10 for ls in values(loadshapes) if ls["name"] != "3")
        @test all(haskey.([loadshapes["$i"] for i in [3, 8, 9]], "qmult"))
        @test all(haskey.([loadshapes["$i"] for i in [4, 6, 8]], "hour"))
    end

    @testset "arrays from files" begin
        dss = PMD.parse_dss("../test/data/opendss/test2_master.dss")

        @test isa(dss["load"]["ld3"]["yearly"], Vector)
        @test isa(dss["load"]["ld2"]["daily"], Vector)

        @test length(dss["load"]["ld3"]["yearly"]) == 10
        @test length(dss["load"]["ld2"]["daily"]) == 10
    end

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

    # TODO fix, do we support these previously erroring cases now?
    # @testset "opendss parse load model errors" begin
    #     dss = PMD.parse_dss("../test/data/opendss/loadparser_error.dss")
    #     for (name, load) in dss["load"]
    #         _dss = deepcopy(dss)
    #         _dss["load"] = Dict{String,Any}(name => load)

    #        @test_throws(TESTLOG, AssertionError, PMD.parse_opendss(_dss))
    #     end
    # end

    @testset "opendss parse load model warnings" begin
        for model in [3, 4, 7, 8]
           dss = PMD.parse_dss("../test/data/opendss/loadparser_warn_model.dss")
           dss["load"] = Dict{String,Any}((n,l) for (n,l) in dss["load"] if l["name"]=="d1phm$model")
           Memento.setlevel!(TESTLOG, "info")
           @test_warn(TESTLOG, ": load model $model not supported. Treating as model 1.", PMD.parse_opendss(dss))
           Memento.setlevel!(TESTLOG, "error")
        end
    end


    @testset "opendss parse generic warnings and errors" begin
        Memento.setlevel!(TESTLOG, "info")

        @test_throws(TESTLOG, ErrorException,
                   PMD.parse_file("../test/data/opendss/test_simple2.dss"))

        @test_warn(TESTLOG, "Command \"solve\" on line 69 in \"test2_master.dss\" is not supported, skipping.",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Command \"show\" on line 71 in \"test2_master.dss\" is not supported, skipping.",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "reactors as constant impedance elements is not yet supported, treating reactor.reactor1 like line",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "line.l1: like=something cannot be found",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Rg,Xg are not fully supported",
                   PMD.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.TestUtils.@test_log(TESTLOG, "info", "Circuit has been reset with the \"clear\" on line 2 in \"test2_master.dss\"",
                               PMD.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.TestUtils.@test_log(TESTLOG, "info", "Redirecting to \"test2_Linecodes.dss\" on line 9 in \"test2_master.dss\"",
                               PMD.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.TestUtils.@test_log(TESTLOG, "info", "Redirecting to \"test2_Loadshape.dss\" on line 10 in \"test2_master.dss\"",
                               PMD.parse_file("../test/data/opendss/test2_master.dss"))

        Memento.setlevel!(TESTLOG, "error")
    end

    eng = PMD.parse_file("../test/data/opendss/test2_master.dss"; data_model="engineering")
    pmd = PMD.parse_file("../test/data/opendss/test2_master.dss")

    len = 0.013516796
    rmatrix=PMD._parse_matrix(Float64, "[1.5000  |0.200000  1.50000  |0.250000  0.25000  2.00000  ]")
    xmatrix=PMD._parse_matrix(Float64, "[1.0000  |0.500000  0.50000  |0.500000  0.50000  1.000000  ]")
    cmatrix = PMD._parse_matrix(Float64, "[8.0000  |-2.00000  9.000000  |-1.75000  -2.50000  8.00000  ]")

    @testset "buscoords automatic parsing" begin
        @test all(haskey(bus, "lon") && haskey(bus, "lat") for bus in values(pmd["bus"]) if "bus_i" in 1:10)
    end

    @testset "opendss parse generic parser verification" begin
        dss = PMD.parse_dss("../test/data/opendss/test2_master.dss")

        @test dss["line"]["l7"]["test_param"] == 100.0

        @test pmd["name"] == "test2"

        @test length(pmd) == 19
        @test length(dss) == 14

        for (key, len) in zip(["bus", "load", "shunt", "branch", "gen", "dcline", "transformer"], [33, 4, 5, 27, 4, 0, 10])
            @test haskey(pmd, key)
            @test length(pmd[key]) == len
        end

        @test all(haskey(dss, key) for key in ["loadshape", "linecode", "buscoords", "options", "filename"])
    end

    # TODO fix, the way we calculate voltage bases changed
    # @testset "opendss parse generic branch values verification" begin
    #     basekv_br5 = pmd["bus"][string(pmd["branch"]["5"]["f_bus"])]["base_kv"]
    #     @test all(isapprox.(pmd["branch"]["5"]["br_r"], rmatrix * len / basekv_br5^2 * pmd["baseMVA"]; atol=1e-6))
    #     @test all(isapprox.(pmd["branch"]["5"]["br_x"], xmatrix * len / basekv_br5^2 * pmd["baseMVA"]; atol=1e-6))
    #     @test all(isapprox.(pmd["branch"]["5"]["b_fr"], diag(basekv_br5^2 / pmd["baseMVA"] * 2.0 * pi * 60.0 * cmatrix * len / 1e9) / 2.0; atol=1e-6))

    #     for i in [6, 3]
    #         basekv_bri = pmd["bus"][string(pmd["branch"]["$i"]["f_bus"])]["base_kv"]
    #         @test all(isapprox.(diag(pmd["branch"]["$i"]["b_fr"]), (3.4 * 2.0 + 1.6) / 3.0 * (basekv_bri^2 / pmd["baseMVA"] * 2.0 * pi * 60.0 / 1e9) / 2.0; atol=1e-6))
    #     end

    #     @test all(isapprox.(pmd["branch"]["9"]["br_r"].*(115/69)^2, diagm(0 => fill(6.3012e-8, 3)); atol=1e-12))
    #     @test all(isapprox.(pmd["branch"]["9"]["br_x"].*(115/69)^2, diagm(0 => fill(6.3012e-7, 3)); atol=1e-12))

    #     for k in ["qd", "pd"]
    #         @test all(isapprox.(pmd["load"]["2"][k], pmd["load"]["2"][k]; atol=1e-12))
    #     end

    #     for k in ["gs", "bs"]
    #         @test all(isapprox.(pmd["shunt"]["1"][k], pmd["shunt"]["3"][k]; atol=1e-12))
    #         @test all(isapprox.(pmd["shunt"]["4"][k], pmd["shunt"]["5"][k]; atol=1e-12))
    #     end

    #     for k in keys(pmd["gen"]["1"])
    #         if !(k in ["gen_bus", "index", "name", "source_id", "connections"])
    #             @test all(isapprox.(pmd["gen"]["3"][k], pmd["gen"]["1"][k]; atol=1e-12))
    #         end
    #     end

    #     for k in keys(pmd["branch"]["11"])
    #         if !(k in ["f_bus", "t_bus", "index", "name", "linecode", "source_id", "t_connections", "f_connections"])
    #             mult = 1.0
    #             if k in ["br_r", "br_x", "g_fr", "g_to", "b_fr", "b_to"]
    #                 # compensation for the different voltage base
    #                 basekv_br5 = pmd["bus"][string(pmd["branch"]["5"]["f_bus"])]["base_kv"]
    #                 basekv_br2 = pmd["bus"][string(pmd["branch"]["2"]["f_bus"])]["base_kv"]
    #                 zmult = (basekv_br5/basekv_br2)^2
    #                 mult = (k in ["br_r", "br_x"]) ? zmult : 1/zmult
    #             end
    #             @test all(isapprox.(pmd["branch"]["5"][k].*mult, pmd["branch"]["2"][k]; atol=1e-12))
    #         end
    #     end
    # end

    @testset "opendss parse length units" begin
        @test eng["line"]["l8"]["length"] == 1000.0 * len
        basekv_br4 = pmd["bus"][string(pmd["branch"]["4"]["f_bus"])]["base_kv"]
        @test all(isapprox.(pmd["branch"]["4"]["br_r"], rmatrix * len / basekv_br4^2 * pmd["baseMVA"]; atol=1e-6))
        @test all(isapprox.(pmd["branch"]["4"]["br_x"], xmatrix * len / basekv_br4^2 * pmd["baseMVA"]; atol=1e-6))
        @test all(isapprox.(pmd["branch"]["4"]["b_fr"], diag(basekv_br4^2 / pmd["baseMVA"] * 2.0 * pi * 60.0 * cmatrix * len / 1e9) / 2.0 / 3; atol=1e-6))
    end

    @testset "opendss parse switch length verify" begin
        @testset "branches with switches" begin
            @test eng["switch"]["_l4"]["length"] == 0.001
            @test !all(get(br, "switch", false) for (_,br) in pmd["branch"] if !startswith(br["name"],"_virtual_branch.switch"))
        end
    end

    @testset "opendss parse transformer parsing verify" begin
        dss_data = PMD.parse_dss("../test/data/opendss/test_transformer_formatting.dss")
        transformer = dss_data["transformer"]["transformer_test"]
        @test transformer["phases"] == 3
        @test transformer["tap"] == PMD._parse_rpn("(0.00625 12 * 1 +)")
        @test transformer["tap_2"] == 1.5
        @test transformer["%loadloss"] == 0.01
        @test transformer["xhl"] == 0.02
        @test transformer["kv_2"] == 12.47
        @test transformer["conn_2"] == "wye"
        @test transformer["tap_3"] == 0.9
        @test transformer["wdg_3"] == 3

        PMD._apply_like!(dss_data["transformer"]["reg4b"], dss_data, "transformer")
        @test dss_data["transformer"]["reg4b"]["%loadloss"] == dss_data["transformer"]["reg4a"]["%loadloss"]

        eng_data = PMD.parse_file("../test/data/opendss/test_transformer_formatting.dss"; data_model="engineering")
        @test all(all(eng_data["transformer"]["$n"]["tm"] .==  tm) for (n, tm) in zip(["transformer_test", "reg4"], [[fill(1.075, 3), fill(1.5, 3), fill(0.9, 3)], [ones(3), ones(3)]]))
    end

    @testset "opendss parse storage" begin
        pmd_storage = PMD.parse_file("../test/data/opendss/case3_balanced_battery.dss")
        for bat in values(pmd_storage["storage"])
            for key in ["energy", "storage_bus", "energy_rating", "charge_rating", "discharge_rating",
                        "charge_efficiency", "discharge_efficiency", "thermal_rating", "qmin", "qmax",
                        "r", "x", "p_loss", "q_loss", "status", "source_id"]
                @test haskey(bat, key)
                if key in ["x", "r", "qmin", "qmax", "thermal_rating"]
                    @test isa(bat[key], Vector)
                end
            end
        end

        @test pmd_storage["storage"]["1"]["source_id"] == "storage.s1"
    end

    @testset "opendss parse verify source_id" begin
        @test pmd["shunt"]["2"]["source_id"] == "capacitor.c1"
        @test pmd["shunt"]["4"]["source_id"] == "reactor.reactor3"

        @test pmd["branch"]["9"]["source_id"] == "line.l1"
        @test pmd["transformer"]["9"]["source_id"] == "_virtual_transformer.transformer.t4.1"  # winding indicated by .1
        @test pmd["branch"]["10"]["source_id"] == "reactor.reactor1"

        @test pmd["gen"]["4"]["source_id"] == "_virtual_gen.vsource.source"
        @test pmd["gen"]["1"]["source_id"] == "generator.g2"

        @test pmd["load"]["1"]["source_id"] == "load.ld1"

        @test all(haskey(component, "source_id") for component_type in PMD._dss_supported_components for component in values(get(pmd, component_type, Dict())) if component_type != "bus")
    end

    @testset "opendss parse verify order of properties on line" begin
        pmd1 = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
        pmd2 = PMD.parse_file("../test/data/opendss/case3_balanced_prop-order.dss")

        delete!(pmd1, "map")
        delete!(pmd2, "map")

        @test pmd1 == pmd2

        dss1 = PMD.parse_dss("../test/data/opendss/case3_balanced.dss")
        dss2 = PMD.parse_dss("../test/data/opendss/case3_balanced_prop-order.dss")

        @test dss1 != dss2
        @test all(a == b for (a, b) in zip(dss2["line"]["ohline"]["prop_order"],["name", "bus1", "bus2", "linecode", "rmatrix", "length"]))
        @test all(a == b for (a, b) in zip(dss2["line"]["quad"]["prop_order"],["name", "bus1", "bus2", "like", "linecode", "length"]))
    end

    @testset "opendss parse verify mvasc3/mvasc1 circuit parse" begin
        dss = PMD.parse_dss("../test/data/opendss/test_simple.dss")
        circuit = PMD._create_vsource("sourcebus", "source"; PMD._to_kwargs(dss["vsource"]["source"])...)

        @test circuit["mvasc1"] == 2100.0
        @test circuit["mvasc3"] == 1900.0
        @test isapprox(circuit["isc3"], 9538.8; atol=1e-1)
        @test isapprox(circuit["isc1"], 10543.0; atol=1e-1)

        dss = PMD.parse_dss("../test/data/opendss/test_simple3.dss")
        circuit = PMD._create_vsource("sourcebus", "source"; PMD._to_kwargs(dss["vsource"]["source"])...)

        @test circuit["mvasc1"] == 2100.0
        @test isapprox(circuit["mvasc3"], 1900.0; atol=1e-1)
        @test circuit["isc3"] == 9538.8
        @test isapprox(circuit["isc1"], 10543.0; atol=1e-1)

        dss = PMD.parse_dss("../test/data/opendss/test_simple4.dss")
        circuit = PMD._create_vsource("sourcebus", "source"; PMD._to_kwargs(dss["vsource"]["source"])...)

        @test isapprox(circuit["mvasc1"], 2091.5; atol=1e-1)
        @test circuit["mvasc3"] == 2000.0
        @test circuit["isc3"] == 10041.0
        @test circuit["isc1"] == 10500.0
    end
end

# @testset "test json parser" begin
#     pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")

#     io = PipeBuffer()
#     PMD.print_file(io, pmd)
#     pmd_json_file = PMD.parse_file(io)

#     @test pmd == pmd_json_file
# end
