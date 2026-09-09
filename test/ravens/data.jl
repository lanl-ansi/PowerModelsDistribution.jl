@info "running misc data handling tests - extended to RAVENS"

@testset "test data handling functions - RAVENS" begin

    @testset "node counting functions" begin
        raven_model_1  = ravens_ut_trans_2w_yy
        @test count_nodes(transform_data_model(raven_model_1)) == 12

        raven_model_1 = ravens_IEEE13_Assets
        @test count_nodes(transform_data_model(raven_model_1)) == 41  # stopped excluding source from node count
    end
end

@testset "test active conductor counting functions - RAVENS" begin
    eng = ravens_ut_trans_2w_yy
    math = transform_data_model(eng)

    n_phases = 3

    math_term_count = count_active_terminals(math)
    math_conn_count = count_active_connections(math)

    @test math_term_count == length(math["bus"]) * n_phases && math_conn_count == (length(math["transformer"]) + length(math["branch"])) * n_phases
end

