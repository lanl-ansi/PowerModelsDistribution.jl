@info "running opendss parser tests"

@testset "test opendss parser" begin
    @testset "bus discovery parsing" begin
        eng = parse_file("../test/data/opendss/test_bus_discovery.dss")

        @test length(eng["bus"]) == 24
        @test all(k in keys(eng["bus"]) for k in [["$i" for i in 1:23]..., "sourcebus"])
    end

    @testset "loadshape parsing" begin
        dss = parse_dss("../test/data/opendss/loadshapes.dss")

        loadshapes = dss.loadshape

        @test isapprox(loadshapes["1"]["interval"], 1.0/60)
        @test all(length(ls["pmult"]) == 10 for ls in values(loadshapes) if ls["name"] != "3")
        @test all(haskey.([loadshapes["$i"] for i in [3, 7, 8]], "qmult"))
        @test all(haskey.([loadshapes["$i"] for i in [4, 6, 8]], "hour"))
    end

    @testset "reverse polish notation" begin
        @test isapprox(PMD._parse_rpn("2 pi * 60 * .001 *"), 2 * pi * 60 * .001; atol=1e-12)
        @test isapprox(PMD._parse_rpn("(14.4 13.8 / sqr 300 *"), (14.4 / 13.8)^2 * 300; atol=1e-12)
        @test isapprox(PMD._parse_rpn("24.9 3 sqrt /"), 24.9 / sqrt(3); atol=1e-12)
        @test all(isapprox.(PMD._parse_array(Float64, "(\"24.9 3 sqrt /\" \"10 2 *\")"), [24.9 / sqrt(3), 10 * 2.0]; atol=1e-12))

        @test PMD._isa_rpn("2 pi * 60 * .001 *")
        @test !PMD._isa_rpn("[ 2 10 ]")

        @test_logs (:warn, "_parse_rpn does not support 'rollup', 'rolldn', or 'swap', leaving as String") match_mode=:any PowerModelsDistribution._parse_rpn("1 2 swap atan2")
        @test_logs (:warn, "' 1 2 + - ' is not valid Reverse Polish Notation, leaving as String") match_mode=:any PowerModelsDistribution._parse_rpn(" 1 2 + - ")
        @test_logs (:warn, "'1 2 3 +' is not valid Reverse Polish Notation, leaving as String") match_mode=:any PowerModelsDistribution._parse_rpn("1 2 3 +")
    end

    @testset "opendss parse load model warnings" begin
        @test_logs (:warn, "dss load model 3 not supported; treating as constant POWER model") (:warn, "dss load model 6 identical to model 1 in current feature set; treating as constant POWER model") (:warn, "dss load model 7 not supported; treating as constant POWER model") (:warn, "dss load model 4 not supported; treating as constant POWER model") match_mode=:any parse_file("../test/data/opendss/loadparser_warn_model.dss")
    end

    @testset "opendss parse spectrum objects" begin
        dss = parse_dss("../test/data/opendss/test2_master.dss")

        for (_, spectrum) in dss["spectrum"]
            @test all(spectrum["harmonic"] .== [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0])
            @test all(spectrum["%mag"] .== [25.0, 0, 50, 0, 75, 0, 100, 0])
            @test all(spectrum["angle"] .== fill(15.0, 8))
        end
    end

    @testset "opendss parse generic warnings and errors" begin
        @test_throws ErrorException parse_file("../test/data/opendss/test_simple2.dss"; data_model=MATHEMATICAL)
        @test_logs (:info, "Command 'solve' on line 70 in 'test2_master.dss' is not supported, skipping.") (:info, "Command 'show' on line 72 in 'test2_master.dss' is not supported, skipping.") (:warn, "reactors as constant impedance elements is not yet supported, treating reactor.reactor1 like line") (:warn, "line.something does not exist, can't apply 'like' on line.l1") (:info, "Circuit has been reset with the 'clear' on line 2 in 'test2_master.dss'") (:info, "Redirecting to 'test2_linecodes.dss' on line 10 in 'test2_Linecodes.dss'") (:info, "Redirecting to 'test2_loadshape.dss' on line 11 in 'test2_Linecodes.dss'") match_mode=:any parse_file("../test/data/opendss/test2_master.dss")
    end

    raw_dss = parse_raw_dss("../test/data/opendss/test2_master.dss")
    dss = parse_dss("../test/data/opendss/test2_master.dss")
    eng = parse_file("../test/data/opendss/test2_master.dss", import_all=true)
    math = parse_file("../test/data/opendss/test2_master.dss"; data_model=MATHEMATICAL, import_all=true)

    @testset "dss . characters in name" begin
        @test haskey(eng["linecode"], "random_.001-test-name")
        @test eng["line"]["test_.002"]["linecode"] == "random_.001-test-name"
    end

    @testset "dss edit command" begin
        @test all(eng["transformer"]["t5"][k] == v for (k,v) in eng["transformer"]["t4"] if !(k in ["name", "bus", "source_id", "rw", "tm_set", "dss", "controls", "tm_nom"]))
        @test all(eng["transformer"]["t5"]["rw"] .== [0.0074, 0.0076])
    end

    @testset "dss assign property in different file" begin
        @test eng["storage"]["s1"]["bus"] == "_b2"
    end

    @testset "subdirectory parsing" begin
        @test haskey(eng["linecode"], "lc7")
    end

    @testset "multiple generation objects on same bus" begin
        @test all(eng["storage"]["s1"]["connections"] .== collect(1:4))
        @test all(eng["solar"]["solar1"]["connections"] .== collect(1:4))
        @test all(eng["bus"]["b1"]["terminals"] .== collect(1:4))
        @test eng["bus"]["b1"]["grounded"] == [4]
    end

    @testset "buscoords automatic parsing" begin
        @test all(haskey(bus, "lon") && haskey(bus, "lat") for bus in values(math["bus"]) if "bus_i" in 1:10)
    end

    @testset "import_all parsing" begin
        @test all(haskey(comp, "dss") && isa(comp["dss"], Dict) for (comp_type, comps) in eng if isa(comps, Dict) for (_,comp) in comps if isa(comp, Dict) && comp_type != "bus" && comp_type != "settings")
        @test all(haskey(comp, "dss") && isa(comp["dss"], Dict) for (comp_type, comps) in math if isa(comps, Dict) for (id,comp) in comps if isa(comp, Dict) && comp_type != "bus" && comp_type != "settings" && comp_type != "map" && !startswith(comp["name"], "_virtual"))
    end

    @testset "opendss parse generic parser verification" begin
        @test math["name"] == "test2"

        @test length(math) == 18
        @test length([p for p in propertynames(raw_dss) if !isempty(getproperty(raw_dss, p))]) == 23

        for (key, len) in zip(["bus", "load", "shunt", "branch", "gen", "transformer", "storage", "switch"], [34, 4, 5, 29, 5, 10, 1, 1])
            @test haskey(math, key)
            @test length(math[key]) == len
        end

        @test all(!isempty(dss[key]) for key in ["loadshape", "linecode", "buscoordinates", "options"])
    end

    @testset "opendss parse matrix and array with mixed delimiters" begin
        @test all(dss["xycurve"]["test_delimiters"]["points"] .== [1.0, 2.0, 3.0, 4.0])
        @test all(dss["linecode"]["test_matrix_syntax"]["rmatrix"] .== [0.1 0 0; 0 0.1 0; 0 0 0.1])
    end

    @testset "dss setbusxy command" begin
        @test dss["buscoordinates"]["testsource"]["x"] == 0.1
        @test dss["buscoordinates"]["testsource"]["y"] == 0.2
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
        @test all(transformer["taps"] .== [1.075, 1.5, 0.9])
        @test all(transformer["kvs"] .== [67.0, 12.47, 12.47])
        @test all(transformer["conns"] .== [DELTA, WYE, WYE])
        @test transformer["%loadloss"] == 0.01
        @test transformer["xhl"] == 0.02

        @test dss_data["transformer"]["reg4b"]["%loadloss"] == dss_data["transformer"]["reg4a"]["%loadloss"]

        eng_data = parse_file("../test/data/opendss/test_transformer_formatting.dss")
        @test all(all(eng_data["transformer"]["$n"]["tm_set"] .==  tm) for (n, tm) in zip(["transformer_test", "reg4"], [[fill(1.075, 3), fill(1.5, 3), fill(0.9, 3)], [ones(3), ones(3)]]))
    end

    @testset "opendss parse line parsing wires - spacing properties" begin
        dss_data = parse_dss("../test/data/opendss/test2_master.dss")

        @test isa(dss_data["line"]["l9"]["wires"], Vector{String}) && all(dss_data["line"]["l9"]["wires"] .== ["wire1", "wire2", "wire1", "wire2"])
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
        @test all(a[1] == b for (a, b) in zip(dss2["line"]["ohline"]["raw_dss"],["name", "bus1", "bus2", "linecode", "rmatrix", "length"]))
        @test all(a[1] == b for (a, b) in zip(dss2["line"]["quad"]["raw_dss"],["name", "like", "bus1", "bus2", "linecode", "length"]))
    end

    @testset "opendss parse verify mvasc3/mvasc1 circuit parse" begin
        dss_data = parse_dss("../test/data/opendss/test_simple.dss")
        circuit = dss_data.vsource["source"]

        @test circuit["mvasc1"] == 2100.0
        @test circuit["mvasc3"] == 1900.0
        @test isapprox(circuit["isc3"], 9538.8; atol=1e-1)
        @test isapprox(circuit["isc1"], 10543.0; atol=1e-1)

        dss_data = parse_dss("../test/data/opendss/test_simple3.dss")
        circuit = dss_data.vsource["source"]

        @test circuit["mvasc1"] == 2100.0
        @test isapprox(circuit["mvasc3"], 1900.0; atol=1e-1)
        @test circuit["isc3"] == 9538.8
        @test isapprox(circuit["isc1"], 10543.0; atol=1e-1)

        dss_data = parse_dss("../test/data/opendss/test_simple4.dss")
        circuit = dss_data.vsource["source"]

        @test isapprox(circuit["mvasc1"], 2091.5; atol=1e-1)
        @test circuit["mvasc3"] == 2000.0
        @test circuit["isc3"] == 10041.0
        @test circuit["isc1"] == 10500.0
    end

    @testset "opendss capcontrol parse" begin
        defaults = filter(x->x.first!="raw_dss",PowerModelsDistribution._convert_model_to_dict(dss["capcontrol"]["c1_ctrl"]))

        @test isempty(PowerModelsDistribution._check_equal(defaults, Dict{String,Any}(
            "name" => "c1_ctrl",
            "element" => "line.l2",
            "capacitor" => "c1",
            "type" => CAP_REACTIVE_POWER,
            "ctphase" => 1,
            "ctratio" => 1.0,
            "deadtime" => 300.0,
            "delay" => 100.0,
            "delayoff" => 100.0,
            "eventlog" => true,
            "offsetting" => -225.0,
            "onsetting" => 150.0,
            "ptphase" => 1,
            "ptratio" => 1.0,
            "terminal" => 1,
            "vbus" => "",
            "vmax" => 7740.0,
            "vmin" => 7110.0,
            "voltoverride" => true,
            "pctminkvar" => 50.0,
            "enabled" => ENABLED,
            "like" => "",
        )))
    end

    @testset "opendss regcontrol parse" begin
        defaults = filter(x->x.first!="raw_dss",PowerModelsDistribution._convert_model_to_dict(dss["regcontrol"]["t1"]))

        @test isempty(PowerModelsDistribution._check_equal(defaults, Dict{String,Any}(
            "name" => "t1",
            "transformer" => "t1",
            "winding" => 2,
            "vreg" => 122.0,
            "band" => 2.0,
            "delay" => 15.0,
            "ptratio" => 20.0,
            "ctprim" => 700.0,
            "r" => 3.0,
            "x" => 9.0,
            "ptphase" => 1,
            "tapwinding" => 2,
            "bus" => "",
            "debugtrace" => false,
            "eventlog" => true,
            "inversetime" => false,
            "maxtapchange" => 16,
            "revband" => 3.0,
            "revdelay" => 60.0,
            "reversible" => false,
            "revneutral" => false,
            "revr" => 0.0,
            "revthreshold" => 100.0,
            "revvreg" => 120.0,
            "revx" => 0.0,
            "tapdelay" => 2.0,
            "tapnum" => 0,
            "vlimit" => 0.0,
            "rev_z" => [0.0, 0.0],
            "ldc_z" => [0.0, 0.0],
            "cogen" => false,
            "remoteptratio" => 60.0,
            "enabled" => ENABLED,
            "like" => "",
        )))
    end

    @testset "tabulation parse" begin
        eng = parse_file("../test/data/opendss/case_tabulation.dss")

        number_buses = 2
        @test length(eng["bus"]) == number_buses
    end
end

@testset "test different regcontrol configurations" begin
    eng = parse_file("../test/data/opendss/IEEE13_test_controls.dss")

    @test all(isequal(eng["transformer"]["reg1"]["controls"]["vreg"][2], [0.0, 118.0, 0.0]))
    @test all(isequal(eng["transformer"]["reg1"]["controls"]["ptratio"][2], [0.0, 22.0, 0.0]))
    @test all(isequal(eng["transformer"]["reg1"]["controls"]["band"][1], [0.0, 0.0, 4.0]))
    @test all(isequal(eng["transformer"]["reg1"]["controls"]["ctprim"][1], [0.0, 0.0, 695.0]))
    @test all(isequal(eng["transformer"]["sub"]["controls"]["r"][1], [0.0, 0.0, 3.0]))
    @test all(isequal(eng["transformer"]["sub"]["controls"]["x"][1], [0.0, 0.0, 9.0]))
end

@testset "test different capcontrol configurations" begin
    eng = parse_file("../test/data/opendss/IEEE13_CapControl.dss")

    @test all(isequal(eng["shunt"]["c1"]["controls"]["type"], CAP_REACTIVE_POWER))
    @test all(isequal(eng["shunt"]["c1"]["controls"]["ptratio"], 1.0))
    @test all(isequal(eng["shunt"]["c2"]["controls"]["type"], [CAP_DISABLED, CAP_VOLTAGE, CAP_VOLTAGE]))
    @test all(isequal(eng["shunt"]["c2"]["controls"]["terminal"], [0, 1, 2]))
    @test all(isequal(eng["shunt"]["c2"]["controls"]["element"], "line.650632"))
    @test all(isequal(eng["shunt"]["c2"]["controls"]["voltoverride"], [false, false, false]))
end

@testset "test json parser" begin
    eng = parse_file("../test/data/opendss/case3_balanced.dss")

    io = PipeBuffer()
    print_file(io, eng)
    eng_json_file = parse_file(io)

    @test eng == eng_json_file
end
