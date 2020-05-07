@info "running data model creation and conversion tests"

@testset "engineering data model" begin
    @testset "helper functions for building engineering data model" begin
        eng = Model()

        add_bus!(eng, "sourcebus"; terminals=[1,2,3,4], grounded=[4])
        add_bus!(eng, "primary"; terminals=[1,2,3])
        add_bus!(eng, "loadbus"; terminals=[1,2,3,4], grounded=[4])

        add_voltage_source!(eng, "source", "sourcebus", [1,2,3,4]; vm=[1, 1, 1])

        add_line!(eng, "trunk", "sourcebus", "primary", [1,2,3], [1,2,3])
        add_line!(eng, "primary", "primary", "loadbus", [1,2,3], [1,2,3])

        add_load!(eng, "balanced", "loadbus", [1,2,3,4]; pd_nom=[5, 5, 5], qd_nom=[1, 1, 1])

        add_vbase_default!(eng, "sourcebus", 1)

        result = run_mc_opf(eng, ACPPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0150; atol=1e-4)

        eng2 = deepcopy(eng)

        add_bus!(eng2, "ttbus"; terminals=[1,2,3,4], grounded=[4])

        add_transformer!(eng2, "tx1", "sourcebus", "ttbus", [1,2,3,4], [1,2,3,4])

        add_bus!(eng2, "loadbus2"; terminals=[1,2,3,4], grounded=[4])

        add_switch!(eng2, "breaker", "ttbus", "loadbus2", [1,2,3], [1,2,3])

        add_load!(eng2, "tload", "loadbus2", [1,2,3,4]; pd_nom=[5, 5, 5], qd_nom=[1, 1, 1])

        add_generator!(eng2, "secondary", "loadbus2", [1,2,3,4]; cost_pg_parameters=[0.0, 1.2, 0])

        add_shunt!(eng2, "cap", "loadbus2", [1,2,3,4]; bs=diagm(0=>fill(1, 3)))

        result2 = run_mc_opf(eng2, ACPPowerModel, ipopt_solver)

        @test result2["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result2["objective"], -83.3003; atol=0.2)
    end

    @testset "engineering model transformations" begin
        eng = parse_file("../test/data/opendss/case3_balanced.dss"; transformations=[(apply_voltage_bounds!, "vm_ub"=>Inf)])

        @test all(all(isapprox.(bus["vm_lb"], 0.4 / sqrt(3) * 0.9)) && all(isinf.(bus["vm_ub"])) for (id,bus) in eng["bus"] if id != "sourcebus")

        math = transform_data_model(eng)

        @test all(all(isapprox.(bus["vmin"], 0.9)) for (_,bus) in math["bus"] if bus["name"] != "sourcebus" && !startswith(bus["name"], "_virtual"))

        eng = parse_file("../test/data/opendss/ut_trans_2w_yy.dss"; transformations=[make_lossless!, (apply_voltage_bounds!, "vm_lb"=>0.95, "vm_ub"=>1.05)])

        @test all(all(eng["transformer"]["tx1"][k] .== 0) for k in ["rw", "xsc", "noloadloss", "imag"])

        math = transform_data_model(eng)

        @test length(math["bus"]) == 5
        @test all(all(isapprox.(bus["vmin"], 0.95)) && all(isapprox.(bus["vmax"], 1.05)) for (_,bus) in math["bus"] if bus["name"] != "sourcebus" && !startswith(bus["name"], "_virtual"))
    end
end
