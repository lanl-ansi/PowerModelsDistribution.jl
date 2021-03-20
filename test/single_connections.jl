@info "running single connections (sc) tests"

@testset "test sc" begin

    case3_unbalanced = parse_file("../test/data/opendss/case3_unbalanced.dss", data_model = MATHEMATICAL)
   
    # remove line charging
    for (b, branch) in case3_unbalanced["branch"] branch["b_to"] = fill(0.0, (3,3)) end
    for (b, branch) in case3_unbalanced["branch"] branch["b_fr"] = fill(0.0, (3,3)) end
    for l in ["2", "3"] delete!(case3_unbalanced["load"], l) end #remaining load "1" is connected to phase 2
    case3_unbalanced["branch"]["1"]["br_x"] = Array(Diagonal(case3_unbalanced["branch"]["1"]["br_x"]))
    case3_unbalanced["branch"]["1"]["br_r"] = Array(Diagonal(case3_unbalanced["branch"]["1"]["br_r"]))

    case3_reduced = deepcopy(case3_unbalanced)
    case3_reduced["bus"]["3"]["terminals"] = [2]
    case3_reduced["bus"]["3"]["grounded"] = Bool[0]
    case3_reduced["bus"]["3"]["vmin"] = [0.0]
    case3_reduced["bus"]["3"]["vmax"] = [Inf]
    case3_reduced["branch"]["1"]["t_connections"], case3_reduced["branch"]["1"]["f_connections"] = [2], [2]
    case3_reduced["branch"]["1"]["br_r"] = [case3_reduced["branch"]["1"]["br_r"][2,2]]
    case3_reduced["branch"]["1"]["br_x"] = [case3_reduced["branch"]["1"]["br_x"][2,2]]
    case3_reduced["branch"]["1"]["g_to"] = [case3_reduced["branch"]["1"]["g_to"][1]]
    case3_reduced["branch"]["1"]["g_fr"] = [case3_reduced["branch"]["1"]["g_fr"][1]]
    case3_reduced["branch"]["1"]["b_to"] = [case3_reduced["branch"]["1"]["b_to"][1]]
    case3_reduced["branch"]["1"]["b_fr"] = [case3_reduced["branch"]["1"]["b_fr"][1]]

    @testset "sc acp pf" begin
        sol = solve_mc_pf(case3_unbalanced, ACPPowerModel, ipopt_solver)
        sol_reduced = solve_mc_pf(case3_reduced, ACPPowerModel, ipopt_solver)
        for b in ["1", "2"]
            @test all(isapprox.(sol["solution"]["bus"][b]["vm"], sol_reduced["solution"]["bus"][b]["vm"]))
            @test all(isapprox.(sol["solution"]["bus"][b]["va"], sol_reduced["solution"]["bus"][b]["va"]))
        end
        @test isapprox(sol["solution"]["bus"]["3"]["vm"][2], sol_reduced["solution"]["bus"]["3"]["vm"][1])
        @test isapprox(sol["solution"]["bus"]["3"]["va"][2], sol_reduced["solution"]["bus"]["3"]["va"][1])
    end

    @testset "sc acr pf" begin
        sol = solve_mc_pf(case3_unbalanced, ACRPowerModel, ipopt_solver)
        sol_reduced = solve_mc_pf(case3_reduced, ACRPowerModel, ipopt_solver)
        for b in ["1", "2"]
            @test all(isapprox.(sol["solution"]["bus"][b]["vr"], sol_reduced["solution"]["bus"][b]["vr"]))
            @test all(isapprox.(sol["solution"]["bus"][b]["vi"], sol_reduced["solution"]["bus"][b]["vi"]))
        end
        @test isapprox(sol["solution"]["bus"]["3"]["vr"][2], sol_reduced["solution"]["bus"]["3"]["vr"][1])
        @test isapprox(sol["solution"]["bus"]["3"]["vi"][2], sol_reduced["solution"]["bus"]["3"]["vi"][1])
    end

    @testset "sc ivr pf" begin
        sol = solve_mc_pf(case3_unbalanced, IVRPowerModel, ipopt_solver)
        sol_reduced = solve_mc_pf(case3_reduced, IVRPowerModel, ipopt_solver)
        for b in ["1", "2"]
            @test all(isapprox.(sol["solution"]["bus"][b]["vr"], sol_reduced["solution"]["bus"][b]["vr"]))
            @test all(isapprox.(sol["solution"]["bus"][b]["vi"], sol_reduced["solution"]["bus"][b]["vi"]))
        end
        @test isapprox(sol["solution"]["bus"]["3"]["vr"][2], sol_reduced["solution"]["bus"]["3"]["vr"][1])
        @test isapprox(sol["solution"]["bus"]["3"]["vi"][2], sol_reduced["solution"]["bus"]["3"]["vi"][1])
    end

    @testset "sc LinDist3Flow pf" begin
        sol = solve_mc_pf(case3_unbalanced, LPUBFDiagPowerModel, ipopt_solver)
        sol_reduced = solve_mc_pf(case3_reduced, LPUBFDiagPowerModel, ipopt_solver)
        for b in ["1", "2"]
            @test all(isapprox.(sol["solution"]["bus"][b]["w"], sol_reduced["solution"]["bus"][b]["w"]))
        end
        @test isapprox(sol["solution"]["bus"]["3"]["w"][2], sol_reduced["solution"]["bus"]["3"]["w"][1])
    end
end