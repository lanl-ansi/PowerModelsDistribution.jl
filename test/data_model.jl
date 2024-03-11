@info "running data model creation and conversion tests"

@testset "engineering data model" begin
    @testset "helper functions for building engineering data model" begin
        eng = Model()

        add_bus!(eng, "sourcebus"; terminals=[1,2,3,4], grounded=[4])
        add_bus!(eng, "primary"; terminals=[1,2,3])
        add_bus!(eng, "loadbus"; terminals=[1,2,3,4], grounded=[4])

        add_voltage_source!(eng, "source", "sourcebus", [1,2,3,4]; vm=[1, 1, 1, 0])

        add_linecode!(eng, "default", LinearAlgebra.diagm(0=>fill(0.01, 3)), LinearAlgebra.diagm(0=>fill(0.2, 3)))

        add_line!(eng, "trunk", "sourcebus", "primary", [1,2,3], [1,2,3]; linecode="default")
        add_line!(eng, "primary", "primary", "loadbus", [1,2,3], [1,2,3]; linecode="default")

        add_load!(eng, "balanced", "loadbus", [1,2,3,4]; pd_nom=[5, 5, 5], qd_nom=[1, 1, 1])

        add_vbase_default!(eng, "sourcebus", 1)

        result = solve_mc_opf(eng, ACPUPowerModel, ipopt_solver)

        @test result["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result["objective"], 0.0150; atol=1e-4)

        eng2 = deepcopy(eng)

        add_bus!(eng2, "ttbus"; terminals=[1,2,3,4], grounded=[4])

        add_transformer!(eng2, "tx1", "sourcebus", "ttbus", [1,2,3,4], [1,2,3,4])

        add_bus!(eng2, "loadbus2"; terminals=[1,2,3,4], grounded=[4])

        add_switch!(eng2, "breaker", "ttbus", "loadbus2", [1,2,3], [1,2,3])

        add_load!(eng2, "tload", "loadbus2", [1,2,3,4]; pd_nom=[5, 5, 5], qd_nom=[1, 1, 1])

        add_generator!(eng2, "secondary", "loadbus2", [1,2,3,4]; cost_pg_parameters=[0.0, 1.2, 0], pg_ub=[2.5, 2.5, 2.5, 0.0], pg_lb=zeros(4), qg_lb=[-2.5, -2.5, -2.5, 0], qg_ub=[2.5, 2.5, 2.5, 0])

        add_shunt!(eng2, "cap", "loadbus2", [1,2,3,4]; bs=LinearAlgebra.diagm(0=>fill(0.1, 3)))

        result2 = solve_mc_opf(eng2, ACRUPowerModel, ipopt_solver)

        @test result2["termination_status"] == LOCALLY_SOLVED
        @test isapprox(result2["objective"], 0.03; atol=0.2)
    end

    @testset "engineering model transformations" begin
        eng = parse_file("../test/data/opendss/case3_balanced.dss"; transformations=[(apply_voltage_bounds!, "vm_ub"=>Inf)])

        @test all(all(isapprox.(bus["vm_lb"][.!any.(bus["grounded"] .== t for t in bus["terminals"])], 0.4 / sqrt(3) * 0.9)) && all(isinf.(bus["vm_ub"])) for (id,bus) in eng["bus"] if id != "sourcebus")

        math = transform_data_model(eng)

        @test all(all(isapprox.(bus["vmin"], 0.9)) for (_,bus) in math["bus"] if bus["name"] != "sourcebus" && !startswith(bus["name"], "_virtual"))

        eng = parse_file("../test/data/opendss/ut_trans_2w_yy.dss"; transformations=[(make_lossless!, "exclude"=>String["line", "linecode"]), (apply_voltage_bounds!, "vm_lb"=>0.95, "vm_ub"=>1.05)])

        @test all(all(eng["transformer"]["tx1"][k] .== 0) for k in ["rw", "xsc", "noloadloss", "cmag"])

        math = transform_data_model(eng)

        @test length(math["bus"]) == 5
        @test all(all(isapprox.(bus["vmin"], 0.95)) && all(isapprox.(bus["vmax"], 1.05)) for (_,bus) in math["bus"] if bus["name"] != "sourcebus" && !startswith(bus["name"], "_virtual"))

        eng = parse_file("../test/data/opendss/case3_balanced.dss"; transformations=[remove_all_bounds!])
        @test !all(haskey(line, "cm_ub") for line in values(eng["line"]))
    end

    @testset "jump model from engineering data model" begin
        eng = parse_file("../test/data/opendss/case3_balanced.dss")

        pm_eng = instantiate_mc_model(eng, ACPUPowerModel, build_mc_opf)

        math = transform_data_model(eng)

        pm_math = instantiate_mc_model(math, ACPUPowerModel, build_mc_opf)

        @test sprint(print, pm_eng) == sprint(print, pm_math)
    end

    @testset "test manual 1-phase model creation" begin
        data = parse_file("../test/data/opendss/case3_unbalanced.dss"; data_model = MATHEMATICAL)
        delete!(data["load"], "2")
        delete!(data["load"], "3")
        data["bus"]["3"]["terminals"] = Vector{Int}([2])
        data["bus"]["3"]["grounded"] = Vector{Bool}([0])

        for (k,v) in data["branch"]["1"]
            if isa(v, Matrix)
                data["branch"]["1"][k] = reshape([v[2,2]], 1, 1)
            elseif isa(v, Vector)
                data["branch"]["1"][k] = [v[2]]
            end
        end

        result_pf = solve_mc_pf(data, ACRUPowerModel, ipopt_solver)
        @test result_pf["objective"] == 0.0

        result_opf = solve_mc_opf(data, ACRUPowerModel, ipopt_solver)
        @test isapprox(result_opf["objective"], 0.006197; atol=1e-3)
    end

    @testset "test reduce_line_series" begin
        eng = parse_file("../test/data/opendss/line_series.dss")
        engn = deepcopy(eng)
        reduce_line_series!(engn)

        r = solve_mc_opf(eng, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])
        rn = solve_mc_opf(engn, ACRUPowerModel, ipopt_solver; solution_processors=[sol_data_model!])

        for bus_id in [line[bus_end] for line in values(engn["line"]) for bus_end in ["f_bus", "t_bus"] if line["status"]!=DISABLED]
            @test all(isapprox.(rn["solution"]["bus"][bus_id]["vm"], r["solution"]["bus"][bus_id]["vm"]; atol=1e-4))
        end
    end
end
