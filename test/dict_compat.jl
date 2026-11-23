@info "running legacy Dict interface compatibility tests"

@testset "legacy Dict interface compatibility" begin
    @testset "dss2eng_extensions receive a Dict" begin
        received_type = Ref{Any}(nothing)
        function _test_dss2eng_extension!(data_eng::Dict{String,<:Any}, data_dss)
            received_type[] = typeof(data_eng)
            data_eng["test_extension_key"] = true
        end

        eng = parse_file("../test/data/opendss/case3_balanced.dss"; dss2eng_extensions=Function[_test_dss2eng_extension!])

        @test received_type[] <: Dict
        @test eng["test_extension_key"]
    end

    @testset "calc_voltage_bases on ENGINEERING Dict" begin
        eng = deepcopy(case3_unbalanced)
        eng_dict = PMD._convert_model_to_dict(eng)

        bus_vbases, _ = calc_voltage_bases(eng_dict, eng_dict["settings"]["vbases_default"])
        bus_vbases_model, _ = calc_voltage_bases(eng, eng["settings"]["vbases_default"])

        @test !any(ismissing, values(bus_vbases))
        @test bus_vbases == bus_vbases_model

        apply_voltage_bounds!(eng_dict; vm_lb=0.9, vm_ub=1.1)
        source_buses = Set(vs["bus"] for vs in values(eng_dict["voltage_source"]))
        @test all(haskey(bus, "vm_lb") && haskey(bus, "vm_ub") for (id, bus) in eng_dict["bus"] if !(id in source_buses))
    end

    @testset "_convert_model_to_dict on Dict-backed models" begin
        eng = deepcopy(case3_unbalanced)
        eng_dict = PMD._convert_model_to_dict(eng)

        @test !haskey(eng_dict, "data")
        @test iseng(eng_dict)
        @test eng_dict["settings"] == eng["settings"]

        mn_eng_dict = PMD._convert_model_to_dict(make_multinetwork(eng))
        @test ismultinetwork(mn_eng_dict)
        @test all(nw isa Dict{String,Any} for nw in values(mn_eng_dict["nw"]))
        @test all(haskey(nw, "settings") for nw in values(mn_eng_dict["nw"]))
    end

    @testset "model equality" begin
        @test deepcopy(case3_unbalanced) == deepcopy(case3_unbalanced)
        @test case3_unbalanced != case3_balanced
    end

    @testset "solve_mc_model with Dict input returns Dict solution" begin
        eng_dict = PMD._convert_model_to_dict(deepcopy(case3_unbalanced))
        result = solve_mc_opf(eng_dict, ACPUPowerModel, ipopt_solver)
        result_model = solve_mc_opf(deepcopy(case3_unbalanced), ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test result["solution"] isa Dict{String,Any}
        @test isapprox(result["objective"], result_model["objective"]; atol=1e-6)

        math_dict = PMD._convert_model_to_dict(transform_data_model(deepcopy(case3_unbalanced)))
        result_math = solve_mc_opf(math_dict, ACPUPowerModel, ipopt_solver)

        @test result_math["termination_status"] == LOCALLY_SOLVED
        @test result_math["solution"] isa Dict{String,Any}
        @test isapprox(result_math["objective"], result_model["objective"]; atol=1e-6)
    end
end
