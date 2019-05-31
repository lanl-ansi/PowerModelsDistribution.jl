@testset "test idempotent units transformations" begin
    @testset "5-bus case" begin
        data = TPPMs.parse_file("../test/data/matlab/case5_i_r_b.m")
        data_base = deepcopy(data)

        PMs.make_mixed_units(data)
        PMs.make_per_unit(data)
        @test InfrastructureModels.compare_dict(data, data_base)
    end
end


@testset "angle wrapper functions" begin
    @testset "wraptopi" begin
        wrappedradians = TPPMs.wraptopi([0, pi/2, pi, 3pi/2, 2pi])
        @test isapprox(wrappedradians, [0, pi/2, -pi, -pi/2, 0]; atol=1e-12)
    end

    @testset "wrapto180" begin
        wrappeddegrees = TPPMs.wrapto180([0, 90, 180, 270, 360])
        @test isapprox(wrappeddegrees, [0, 90, -180, -90, 0]; atol=1e-12)
    end
end

@testset "matrix manipulation functions" begin
    @testset "3x3 matrix" begin
        A = [1 2 3; 4 5 6; 7 8.0 9]
        utrivec = TPPMs.mat2utrivec(A)
        @test isapprox(utrivec, [2, 3, 6])
        ltrivec = TPPMs.mat2ltrivec(A)
        @test isapprox(ltrivec, [4, 7, 8])
        @test isapprox(A, diagm(0 => diag(A)) + TPPMs.vec2utri(utrivec) + TPPMs.vec2ltri(ltrivec))
    end
    @testset "5x5 matrix" begin
        A = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20; 21 22 23 24 25]
        n = size(A, 1)
        utrivec = TPPMs.mat2utrivec(A)
        @test isapprox(length(utrivec), (n^2-n)/2)
        ltrivec = TPPMs.mat2ltrivec(A)
        @test isapprox(A, diagm(0 => diag(A)) + TPPMs.vec2utri(utrivec) + TPPMs.vec2ltri(ltrivec))
    end
    @testset "3x3 Hermitian matrix" begin
        A = [1 2 3; 4 5 6; 7 8.0 9]
        (r, i) = TPPMs.make_hermitian_matrix_variable(diag(A), TPPMs.mat2utrivec(A), TPPMs.mat2utrivec(A))
        @test size(r) == size(A)
        @test issymmetric(r)
        @test ishermitian(r + im*i)
    end

    @testset "5x5 Hermitian matrix" begin
        A = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20; 21 22 23 24 25]
        (r, i) = TPPMs.make_hermitian_matrix_variable(diag(A), TPPMs.mat2utrivec(A), TPPMs.mat2utrivec(A))
        @test size(r) == size(A)
        @test issymmetric(r)
        @test ishermitian(r + im*i)
    end

    @testset "5x5 full matrix" begin
        A = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20; 21 22 23 24 25]
        B = TPPMs.make_full_matrix_variable(diag(A), TPPMs.mat2ltrivec(A), TPPMs.mat2utrivec(A))
        @test isapprox(A, B)
    end
end
