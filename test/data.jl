@info "running misc data handling tests"

@testset "test impedance to admittance" begin
    branch = Dict{String, Any}()
    branch["br_r"] = [1 2;3 4]
    branch["br_x"] = [1 2;3 4]
    g,b  = calc_branch_y(branch)

    @test typeof(g) <: Matrix
    @test isapprox(g, [-1.0 0.5; 0.75 -0.25])
    @test isapprox(b, [1.0 -0.5; -0.75 0.25])

    branch["br_r"] = [1 2 0;3 4 0; 0 0 0]
    branch["br_x"] = [1 2 0;3 4 0; 0 0 0]
    g,b  = calc_branch_y(branch)

    @test typeof(g) <: Matrix
    @test isapprox(g, [-1.0 0.5 0; 0.75 -0.25 0; 0 0 0])
    @test isapprox(b, [1.0 -0.5 0; -0.75 0.25 0; 0 0 0])
end

@testset "test data handling functions" begin

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
        dss  = parse_dss("../test/data/opendss/case5_phase_drop.dss")
        eng  = parse_file("../test/data/opendss/case5_phase_drop.dss")
        math = parse_file("../test/data/opendss/case5_phase_drop.dss"; data_model=MATHEMATICAL)

        @test count_nodes(dss) == 10 # stopped excluding source from node count
        @test count_nodes(dss) == count_nodes(eng)
        @test count_nodes(eng) == count_nodes(math)

        dss = parse_dss("../test/data/opendss/ut_trans_2w_yy.dss")
        @test count_nodes(dss) == 12  # stopped excluding source from node count
    end
end

@testset "test active conductor counting functions" begin
    eng = parse_file("../test/data/opendss/ut_trans_2w_yy.dss")
    math = transform_data_model(eng)

    n_phases = 3

    eng_term_count = count_active_terminals(eng)
    eng_conn_count = count_active_connections(eng)

    @test eng_term_count == length(eng["bus"]) * n_phases && eng_conn_count == (length(eng["transformer"]) + length(eng["line"])) * n_phases

    math_term_count = count_active_terminals(math)
    math_conn_count = count_active_connections(math)

    @test math_term_count == length(math["bus"]) * n_phases && math_conn_count == (length(math["transformer"]) + length(math["branch"])) * n_phases
end
