@info "running misc data handling tests"

@testset "test impedance to admittance" begin
    branch["br_r"] = PowerModels.MultiConductorMatrix([1 2;3 4])
    branch["br_x"] = PowerModels.MultiConductorMatrix([1 2;3 4])
    g,b  = PowerModels.calc_branch_y(branch)

    @test typeof(g) <: PowerModels.MultiConductorMatrix
    @test isapprox(g.values, [-1.0 0.5; 0.75 -0.25])
    @test isapprox(b.values, [1.0 -0.5; -0.75 0.25])

    branch["br_r"] = PowerModels.MultiConductorMatrix([1 2 0;3 4 0; 0 0 0])
    branch["br_x"] = PowerModels.MultiConductorMatrix([1 2 0;3 4 0; 0 0 0])
    g,b  = PowerModels.calc_branch_y(branch)

    @test typeof(g) <: PowerModels.MultiConductorMatrix
    @test isapprox(g.values, [-1.0 0.5 0; 0.75 -0.25 0; 0 0 0])
    @test isapprox(b.values, [1.0 -0.5 0; -0.75 0.25 0; 0 0 0])
end

@testset "test data handling functions" begin
    @testset "test idempotent units transformations - 5-bus case" begin
        data = PMD.parse_file("../test/data/matlab/case5_i_r_b.m")
        data_base = deepcopy(data)

        PMs.make_mixed_units!(data)
        PMs.make_per_unit!(data)
        @test data == data_base
    end

    @testset "angle wrapper functions" begin
        wrappedradians = PMD._wrap_to_pi([0, pi/2, pi, 3pi/2, 2pi])
        @test isapprox(wrappedradians, [0, pi/2, -pi, -pi/2, 0]; atol=1e-12)

        wrappeddegrees = PMD._wrap_to_180([0, 90, 180, 270, 360])
        @test isapprox(wrappeddegrees, [0, 90, -180, -90, 0]; atol=1e-12)
    end

    @testset "matrix manipulation functions" begin
        @testset "3x3 matrix manipulation" begin
            A = [1 2 3; 4 5 6; 7 8.0 9]
            utrivec = PMD._mat2utrivec!(A)
            @test isapprox(utrivec, [2, 3, 6])
            ltrivec = PMD._mat2ltrivec!(A)
            @test isapprox(ltrivec, [4, 7, 8])
            @test isapprox(A, diagm(0 => diag(A)) + PMD._vec2utri!(utrivec) + PMD._vec2ltri!(ltrivec))
        end
        @testset "5x5 matrix manipulation" begin
            A = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20; 21 22 23 24 25]
            n = size(A, 1)
            utrivec = PMD._mat2utrivec!(A)
            @test isapprox(length(utrivec), (n^2-n)/2)
            ltrivec = PMD._mat2ltrivec!(A)
            @test isapprox(A, diagm(0 => diag(A)) + PMD._vec2utri!(utrivec) + PMD._vec2ltri!(ltrivec))
        end
        @testset "3x3 Hermitian matrix manipulation" begin
            A = [1 2 3; 4 5 6; 7 8.0 9]
            (r, i) = PMD._make_hermitian_matrix_variable(diag(A), PMD._mat2utrivec!(A), PMD._mat2utrivec!(A))
            @test size(r) == size(A)
            @test issymmetric(r)
            @test ishermitian(r + im*i)
        end

        @testset "5x5 Hermitian matrix manipulation" begin
            A = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20; 21 22 23 24 25]
            (r, i) = PMD._make_hermitian_matrix_variable(diag(A), PMD._mat2utrivec!(A), PMD._mat2utrivec!(A))
            @test size(r) == size(A)
            @test issymmetric(r)
            @test ishermitian(r + im*i)
        end

        @testset "5x5 full matrix manipulation" begin
            A = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20; 21 22 23 24 25]
            B = PMD._make_full_matrix_variable(diag(A), PMD._mat2ltrivec!(A), PMD._mat2utrivec!(A))
            @test isapprox(A, B)
        end
    end

    @testset "node counting functions" begin
        dss = PMD.parse_dss("../test/data/opendss/case5_phase_drop.dss")
        pmd = PMD.parse_file("../test/data/opendss/case5_phase_drop.dss")
        matlab = PMD.parse_file("../test/data/matlab/case5_i_r_b.m")

        @test count_nodes(dss) == 7
        @test count_nodes(dss) == count_nodes(pmd)
        @test count_nodes(matlab) == 15

        dss = PMD.parse_dss("../test/data/opendss/ut_trans_2w_yy.dss")
        @test count_nodes(dss) == 9
    end
end
