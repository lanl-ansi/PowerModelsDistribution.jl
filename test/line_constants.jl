@testset "line constants" begin
    eng = parse_file("../test/data/opendss/IEEE13_Assets.dss")

    z = [
        0.000211485 9.80031e-5  9.61061e-5;
        9.80031e-5  0.000217924 9.93019e-5;
        9.61061e-5  9.93019e-5  0.000213946
    ] + 1im * [
        0.000653882 0.000299058 0.000241464;
        0.000299058 0.000632756 0.000273149;
        0.000241464 0.000273149 0.000645754
    ]
    y = [
        0 0 0;
        0 0 0;
        0 0 0
    ] + 1im * [
         0.00961082 -0.00291343 -0.00123058;
        -0.00291343  0.0103002  -0.00230096;
        -0.00123058 -0.00230096  0.00939874
    ]

    @test all(isapprox.(eng["line"]["650632"]["rs"] + 1im * eng["line"]["650632"]["xs"], z; atol=1e-4))
    @test all(isapprox.(eng["line"]["650632"]["g_fr"].*2 + 1im * eng["line"]["650632"]["b_fr"].*2, y; atol=1e-4))

    z = [
        0.000491885 0.000198743 0.000177491;
        0.000198743 0.000486187 0.000198743;
        0.000177491 0.000198743 0.000491885
    ] + 1im * [
         0.00027695 1.99955e-5  -9.24851e-6;
         1.99955e-5 0.000250698  1.99955e-5;
        -9.24851e-6 1.99955e-5   0.00027695
    ]
    y = [
        0 0 0;
        0 0 0;
        0 0 0
    ] + 1im * [
        0.238581 0.0 0.0;
        0.0 0.238581 0.0;
        0.0 0.0 0.238581
    ]

    @test all(isapprox.(eng["line"]["692675"]["rs"] + 1im * eng["line"]["692675"]["xs"], z; atol=1e-5))
    @test all(isapprox.(eng["line"]["692675"]["g_fr"].*2 + 1im * eng["line"]["692675"]["b_fr"].*2, y; atol=1e-5))

    z = [0.000802907;;] + 1im * [0.000429279;;]
    y = [0;;] + 1im * [0.166359;;]

    @test all(isapprox.(eng["line"]["684652"]["rs"] + 1im * eng["line"]["684652"]["xs"], z; atol=1e-5))
    @test all(isapprox.(eng["line"]["684652"]["g_fr"].*2 + 1im * eng["line"]["684652"]["b_fr"].*2, y; atol=1e-5))

    r = solve_mc_pf(eng, ACRUPowerModel, ipopt_solver; make_si=false, solution_processors=[sol_data_model!])

    @test all(isapprox.(r["solution"]["bus"]["652"]["vm"], [0.974806]; atol=1e-4))
    @test all(isapprox.(r["solution"]["bus"]["652"]["va"], [-5.35]; atol=1e-0))

    @test all(isapprox.(r["solution"]["bus"]["675"]["vm"], [0.975987, 1.04705, 0.984996]; atol=1e-4))
    @test all(isapprox.(r["solution"]["bus"]["675"]["va"], [-5.52, -121.71, 116.03]; atol=1e-0))
end
