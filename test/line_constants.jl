@testset "line constants" begin
    eng = parse_file("../test/data/opendss/IEEE13_Assets.dss")

    z = [0.000211485 9.80031e-5 9.61061e-5; 9.80031e-5 0.000217924 9.93019e-5; 9.61061e-5 9.93019e-5 0.000213946] + 1im * [0.000653882 0.000299058 0.000241464; 0.000299058 0.000632756 0.000273149; 0.000241464 0.000273149 0.000645754]
    y = [0 0 0; 0 0 0; 0 0 0] + 1im * [0.00961103 -0.00291349 -0.00123061; -0.00291349 0.0103004 -0.00230101; -0.00123061 -0.00230101 0.0093989]

    @test all(isapprox.(eng["line"]["650632"]["rs"] + 1im * eng["line"]["650632"]["xs"], z; atol=1e-4))
    @test all(isapprox.(eng["line"]["650632"]["g_fr"].*2 + 1im * eng["line"]["650632"]["b_fr"].*2, y; atol=1e-4))

    z = [0.000495446  0.000201749  0.000180103; 0.000201749 0.000489627 0.000201749; 0.000180103 0.000201749 0.000495446] + 1im * [-0.000280292 1.5444e-5 -1.4004e-5; 1.5444e-5 -0.000307078 1.5444e-5; -1.4004e-5 1.5444e-5 -0.000280292]
    y = [0 0 0; 0 0 0; 0 0 0] + 1im * [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]

    @test all(isapprox.(eng["line"]["692675"]["rs"] + 1im * eng["line"]["692675"]["xs"], z; atol=1e-4))
    @test all(isapprox.(eng["line"]["692675"]["g_fr"].*2 + 1im * eng["line"]["692675"]["b_fr"].*2, y; atol=1e-4))

    z = [0.000805623042762814;;] + 1im * [-0.00012323041308967396;;]
    y = [0;;] + 1im * [0.0;;]

    @test all(isapprox.(eng["line"]["684652"]["rs"] + 1im * eng["line"]["684652"]["xs"], z; atol=1e-4))
    @test all(isapprox.(eng["line"]["684652"]["g_fr"].*2 + 1im * eng["line"]["684652"]["b_fr"].*2, y; atol=1e-4))

    r = solve_mc_pf(eng, ACRUPowerModel, ipopt_solver; make_si=false, solution_processors=[sol_data_model!])

    @test all(isapprox.(r["solution"]["bus"]["652"]["vm"], [0.977265]; atol=1e-4))
    @test all(isapprox.(r["solution"]["bus"]["652"]["va"], [-5.5489]; atol=1e-2))

    @test all(isapprox.(r["solution"]["bus"]["675"]["vm"], [0.976462, 1.04466, 0.985318]; atol=1e-4))
    @test all(isapprox.(r["solution"]["bus"]["675"]["va"], [-5.55622, -121.977, 116.092]; atol=1e-2))
end
