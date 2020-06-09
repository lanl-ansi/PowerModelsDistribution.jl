@info "running opendss parser tests"

@testset "test opendss parser" begin
    @testset "bus discovery parsing" begin
        eng = parse_file("../test/data/opendss/test_bus_discovery.dss")

        @test length(eng["bus"]) == 24
        @test all(k in keys(eng["bus"]) for k in [["$i" for i in 1:23]..., "sourcebus"])
    end

    @testset "loadshape parsing" begin
        dss = parse_dss("../test/data/opendss/loadshapes.dss")

        loadshapes = Dict{String,Any}()
        for (name, ls) in dss["loadshape"]
            loadshapes[name] = PMD._create_loadshape(name; PMD._to_kwargs(ls)...)
        end

        @test isapprox(loadshapes["1"]["interval"], 1.0/60)
        @test all(length(ls["pmult"]) == 10 for ls in values(loadshapes) if ls["name"] != "3")
        @test all(haskey.([loadshapes["$i"] for i in [3, 8, 9]], "qmult"))
        @test all(haskey.([loadshapes["$i"] for i in [4, 6, 8]], "hour"))
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

    @testset "opendss parse load model warnings" begin
        for model in [3, 4, 7, 8]
           dss = parse_dss("../test/data/opendss/loadparser_warn_model.dss")
           dss["load"] = Dict{String,Any}((n,l) for (n,l) in dss["load"] if n=="d1phm$model")
           Memento.setlevel!(TESTLOG, "info")
           @test_warn(TESTLOG, ": dss load model $model not supported. Treating as constant POWER model", parse_opendss(dss))
           Memento.setlevel!(TESTLOG, "error")
        end
    end


    @testset "opendss parse generic warnings and errors" begin
        Memento.setlevel!(TESTLOG, "info")

        @test_throws(TESTLOG, ErrorException,
                   parse_file("../test/data/opendss/test_simple2.dss"; data_model=MATHEMATICAL))

        @test_warn(TESTLOG, "Command \"solve\" on line 69 in \"test2_master.dss\" is not supported, skipping.",
                   parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Command \"show\" on line 71 in \"test2_master.dss\" is not supported, skipping.",
                   parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "reactors as constant impedance elements is not yet supported, treating reactor.reactor1 like line",
                   parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "line.l1: like=something cannot be found",
                   parse_file("../test/data/opendss/test2_master.dss"))

        @test_warn(TESTLOG, "Rg,Xg are not fully supported",
                   parse_file("../test/data/opendss/test2_master.dss"))

        Memento.TestUtils.@test_log(TESTLOG, "info", "Circuit has been reset with the \"clear\" on line 2 in \"test2_master.dss\"",
                               parse_file("../test/data/opendss/test2_master.dss"))

        Memento.TestUtils.@test_log(TESTLOG, "info", "Redirecting to \"test2_Linecodes.dss\" on line 9 in \"test2_master.dss\"",
                               parse_file("../test/data/opendss/test2_master.dss"))

        Memento.TestUtils.@test_log(TESTLOG, "info", "Redirecting to \"test2_Loadshape.dss\" on line 10 in \"test2_master.dss\"",
                               parse_file("../test/data/opendss/test2_master.dss"))

        Memento.setlevel!(TESTLOG, "error")
    end

    eng = parse_file("../test/data/opendss/test2_master.dss", import_all=true)
    math = parse_file("../test/data/opendss/test2_master.dss"; data_model=MATHEMATICAL, import_all=true)

    @testset "buscoords automatic parsing" begin
        @test all(haskey(bus, "lon") && haskey(bus, "lat") for bus in values(math["bus"]) if "bus_i" in 1:10)
    end

    @testset "import_all parsing" begin
        @test all(haskey(comp, "dss") && isa(comp["dss"], Dict) for (comp_type, comps) in eng if isa(comps, Dict) for (_,comp) in comps if isa(comp, Dict) && comp_type != "bus" && comp_type != "settings")
        @test all(haskey(comp, "dss") && isa(comp["dss"], Dict) for (comp_type, comps) in math if isa(comps, Dict) for (id,comp) in comps if isa(comp, Dict) && comp_type != "bus" && comp_type != "settings" && comp_type != "map" && !startswith(comp["name"], "_virtual"))
    end

    @testset "opendss parse generic parser verification" begin
        dss = parse_dss("../test/data/opendss/test2_master.dss")

        @test dss["line"]["l7"]["test_param"] == 100.0

        @test math["name"] == "test2"

        @test length(math) == 20
        @test length(dss) == 18

        for (key, len) in zip(["bus", "load", "shunt", "branch", "gen", "dcline", "transformer"], [33, 4, 5, 28, 4, 0, 10])
            @test haskey(math, key)
            @test length(math[key]) == len
        end

        @test all(haskey(dss, key) for key in ["loadshape", "linecode", "buscoords", "options", "filename"])
    end

    @testset "opendss parse like" begin
        for k in ["pd_nom", "qd_nom"]
            @test all(isapprox.(eng["load"]["ld2"][k], eng["load"]["ld4"][k]; atol=1e-12))
        end

        for (k,v) in eng["generator"]["g2"]
            if isa(v, Real)
                @test all(isapprox.(v, eng["generator"]["g3"][k]; atol=1e-12))
            end
        end
    end

    @testset "opendss parse length units" begin
        @test eng["line"]["l8"]["length"] == 1000.0 * 0.013516796
    end

    @testset "opendss parse switch length verify" begin
        @testset "branches with switches" begin
            @test eng["switch"]["_l4"]["length"] == 0.001
            @test !all(get(br, "switch", false) for (_,br) in math["branch"] if !startswith(br["name"],"_virtual_branch.switch"))
        end
    end

    @testset "opendss parse transformer parsing verify" begin
        dss_data = parse_dss("../test/data/opendss/test_transformer_formatting.dss")
        transformer = dss_data["transformer"]["transformer_test"]
        @test transformer["phases"] == 3
        @test transformer["tap"] == PMD._parse_rpn("(0.00625 12 * 1 +)")
        @test transformer["tap_2"] == 1.5
        @test transformer["%loadloss"] == 0.01
        @test transformer["xhl"] == 0.02
        @test transformer["kv_2"] == 12.47
        @test transformer["conn_2"] == WYE
        @test transformer["tap_3"] == 0.9
        @test transformer["wdg_3"] == 3

        PMD._apply_like!(dss_data["transformer"]["reg4b"], dss_data, "transformer")
        @test dss_data["transformer"]["reg4b"]["%loadloss"] == dss_data["transformer"]["reg4a"]["%loadloss"]

        eng_data = parse_file("../test/data/opendss/test_transformer_formatting.dss")
        @test all(all(eng_data["transformer"]["$n"]["tm_set"] .==  tm) for (n, tm) in zip(["transformer_test", "reg4"], [[fill(1.075, 3), fill(1.5, 3), fill(0.9, 3)], [ones(3), ones(3)]]))
    end

    @testset "opendss parse line parsing wires - spacing properties" begin
        dss_data = parse_dss("../test/data/opendss/test2_master.dss")

        @test isa(dss_data["line"]["l9"]["wires"], Vector{String}) && all(dss_data["line"]["l9"]["wires"] .== ["wire1", "wire2"])
        @test dss_data["line"]["l9"]["spacing"] == "test_spacing"
        @test haskey(dss_data, "wiredata") && (haskey(dss_data["wiredata"], "wire1") && haskey(dss_data["wiredata"], "wire2"))
        @test haskey(dss_data, "linespacing") && haskey(dss_data["linespacing"], "test_spacing")
    end

    @testset "opendss parse storage" begin
        math_storage = parse_file("../test/data/opendss/case3_balanced_battery.dss"; data_model=MATHEMATICAL)
        for bat in values(math_storage["storage"])
            for key in ["energy", "storage_bus", "energy_rating", "charge_rating", "discharge_rating",
                        "charge_efficiency", "discharge_efficiency", "thermal_rating", "qmin", "qmax",
                        "r", "x", "p_loss", "q_loss", "status", "source_id"]
                @test haskey(bat, key)
                if key in ["x", "r", "qmin", "qmax", "thermal_rating"]
                    @test isa(bat[key], Vector)
                end
            end
        end

        @test math_storage["storage"]["1"]["source_id"] == "storage.s1"
    end

    @testset "opendss parse verify source_id" begin
        @test all(haskey(component, "source_id") for component_type in PMD._dss_supported_components for component in values(get(math, component_type, Dict())) if component_type != "bus")
    end

    @testset "opendss parse verify order of properties on line" begin
        math1 = parse_file("../test/data/opendss/case3_balanced.dss"; data_model=MATHEMATICAL)
        math2 = parse_file("../test/data/opendss/case3_balanced_prop-order.dss"; data_model=MATHEMATICAL)

        delete!(math1, "map")
        delete!(math2, "map")

        @test math1 == math2

        dss1 = parse_dss("../test/data/opendss/case3_balanced.dss")
        dss2 = parse_dss("../test/data/opendss/case3_balanced_prop-order.dss")

        @test dss1 != dss2
        @test all(a == b for (a, b) in zip(dss2["line"]["ohline"]["prop_order"],["name", "bus1", "bus2", "linecode", "rmatrix", "length"]))
        @test all(a == b for (a, b) in zip(dss2["line"]["quad"]["prop_order"],["name", "like", "bus1", "bus2", "linecode", "length"]))
    end

    @testset "opendss parse verify mvasc3/mvasc1 circuit parse" begin
        dss = parse_dss("../test/data/opendss/test_simple.dss")
        circuit = PMD._create_vsource("source"; PMD._to_kwargs(dss["vsource"]["source"])...)

        @test circuit["mvasc1"] == 2100.0
        @test circuit["mvasc3"] == 1900.0
        @test isapprox(circuit["isc3"], 9538.8; atol=1e-1)
        @test isapprox(circuit["isc1"], 10543.0; atol=1e-1)

        dss = parse_dss("../test/data/opendss/test_simple3.dss")
        circuit = PMD._create_vsource("source"; PMD._to_kwargs(dss["vsource"]["source"])...)

        @test circuit["mvasc1"] == 2100.0
        @test isapprox(circuit["mvasc3"], 1900.0; atol=1e-1)
        @test circuit["isc3"] == 9538.8
        @test isapprox(circuit["isc1"], 10543.0; atol=1e-1)

        dss = parse_dss("../test/data/opendss/test_simple4.dss")
        circuit = PMD._create_vsource("source"; PMD._to_kwargs(dss["vsource"]["source"])...)

        @test isapprox(circuit["mvasc1"], 2091.5; atol=1e-1)
        @test circuit["mvasc3"] == 2000.0
        @test circuit["isc3"] == 10041.0
        @test circuit["isc1"] == 10500.0
    end
end

@testset "test json parser" begin
    eng = parse_file("../test/data/opendss/case3_balanced.dss")

    io = PipeBuffer()
    print_file(io, eng)
    eng_json_file = parse_file(io)

    @test eng == eng_json_file
end
