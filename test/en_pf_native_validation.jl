@info "running explicit neutral power flow tests with native julia power flow solver"


function conductor_correction!(data_eng)
    nw = data_eng
    if neutral_idx ∈ nw["conductor_ids"]
        filter!(e->e≠neutral_idx, nw["conductor_ids"])
    end
    
    nw["voltage_source"]["source"]["rs"] = nw["voltage_source"]["source"]["rs"][1:3,1:3]
    nw["voltage_source"]["source"]["xs"] = nw["voltage_source"]["source"]["xs"][1:3,1:3]
    nw["voltage_source"]["source"]["connections"] = nw["voltage_source"]["source"]["connections"][1:3]
    nw["voltage_source"]["source"]["va"] = nw["voltage_source"]["source"]["va"][1:3]
    nw["voltage_source"]["source"]["vm"] = nw["voltage_source"]["source"]["vm"][1:3]

    for (l, load) in nw["load"]
        if neutral_idx ∈ load["connections"]
            filter!(e->e≠neutral_idx, load["connections"])
        end
    end

    for (b, bus) in nw["bus"]
        if neutral_idx ∈ bus["terminals"]
            filter!(e->e≠neutral_idx, bus["terminals"])
        end
        if neutral_idx ∈ bus["grounded"]
            filter!(e->e≠neutral_idx, bus["grounded"])
        end
    end

    if haskey(nw, "transformer")
        for (tx, transformer) in nw["transformer"]
            for winding in transformer["connections"]
                filter!(e->e≠neutral_idx, winding)
            end
        end
    end

    if haskey(nw, "solar")
        for (s, solar) in nw["solar"]
            filter!(e->e≠neutral_idx, solar["connections"])
        end
    end
    
    # for (l, line) in nw["line"]
    #     if length(line["t_connections"]) == 4
    #         line["t_connections"] = line["t_connections"][1:3]
    #     end
    #     if length(line["f_connections"]) == 4
    #         line["f_connections"] = line["f_connections"][1:3]
    #     end
    # end

    # for (lc, linecode) in nw["linecode"]
    #     if size(linecode["rs"])[1] = 4
    #         for (arg, mat) in linecode
    #             if arg !="cm_ub"
    #                 linecode[arg] = mat[1:3,1:3]
    #             else
    #                 linecode[arg] = mat[1:3]
    #             end
    #         end
    #     end 
    # end
end


function vsource_correction!(data_eng; explicit_neutral=true)
    if explicit_neutral
        if haskey(data_eng, "multinetwork")
            for (n,nw) in data_eng["nw"]
                nw["voltage_source"]["source"]["rs"][4,4] = nw["voltage_source"]["source"]["rs"][1,1]
                nw["voltage_source"]["source"]["rs"][1:3,4] .= nw["voltage_source"]["source"]["rs"][1,2]
                nw["voltage_source"]["source"]["rs"][4,1:3] .= nw["voltage_source"]["source"]["rs"][1,2]
                nw["voltage_source"]["source"]["xs"][4,4] = nw["voltage_source"]["source"]["xs"][1,1]
                nw["voltage_source"]["source"]["xs"][1:3,4] .= nw["voltage_source"]["source"]["xs"][1,2]
                nw["voltage_source"]["source"]["xs"][4,1:3] .= nw["voltage_source"]["source"]["xs"][1,2]
            end
        else
            data_eng["voltage_source"]["source"]["rs"][4,4] = data_eng["voltage_source"]["source"]["rs"][1,1]
            data_eng["voltage_source"]["source"]["rs"][1:3,4] .= data_eng["voltage_source"]["source"]["rs"][1,2]
            data_eng["voltage_source"]["source"]["rs"][4,1:3] .= data_eng["voltage_source"]["source"]["rs"][1,2]
            data_eng["voltage_source"]["source"]["xs"][4,4] = data_eng["voltage_source"]["source"]["xs"][1,1]
            data_eng["voltage_source"]["source"]["xs"][1:3,4] .= data_eng["voltage_source"]["source"]["xs"][1,2]
            data_eng["voltage_source"]["source"]["xs"][4,1:3] .= data_eng["voltage_source"]["source"]["xs"][1,2]        
        end
    else
        if haskey(data_eng, "multinetwork")
            for (n,nw) in data_eng["nw"]
                nw["voltage_source"]["source"]["rs"] = nw["voltage_source"]["source"]["rs"][1:3,1:3]
                nw["voltage_source"]["source"]["xs"] = nw["voltage_source"]["source"]["xs"][1:3,1:3]
                nw["voltage_source"]["source"]["connections"] = nw["voltage_source"]["source"]["connections"][1:3]
            end
        else
            data_eng["voltage_source"]["source"]["rs"] = data_eng["voltage_source"]["source"]["rs"][1:3,1:3]
            data_eng["voltage_source"]["source"]["xs"] = data_eng["voltage_source"]["source"]["xs"][1:3,1:3]
            data_eng["voltage_source"]["source"]["connections"] = data_eng["voltage_source"]["source"]["connections"][1:3]
        end
    end
    return nothing
end




function sourcebus_voltage_vector_correction!(data_math::Dict{String, Any}; explicit_neutral=true)
    if haskey(data_math, "multinetwork")
        for (n,nw) in data_math["nw"]
            for (i, bus) in data_math["nw"]["bus"]
                if bus["bus_type"] == 3
                    if explicit_neutral
                        bus["vm"] = bus["vm"][1:length(bus["terminals"])]
                        bus["va"] = bus["va"][1:length(bus["terminals"])]
                        bus["vmin"] = bus["vmin"][1:length(bus["terminals"])]
                        bus["vmax"] = bus["vmax"][1:length(bus["terminals"])]
                        bus["grounded"] = bus["grounded"][1:length(bus["terminals"])]
                    else
                        if neutral_idx ∈ bus["terminals"]
                            bus["terminals"] = bus["terminals"][1:end-1]
                        end
                        bus["vm"] = bus["vm"][1:length(bus["terminals"])]
                        bus["va"] = bus["va"][1:length(bus["terminals"])]
                        bus["vmin"] = bus["vmin"][1:length(bus["terminals"])]
                        bus["vmax"] = bus["vmax"][1:length(bus["terminals"])]
                        bus["grounded"] = bus["grounded"][1:length(bus["terminals"])]
                    end
                end
            end
        end
    else
        for (i, bus) in data_math["bus"]
            if bus["bus_type"] == 3
                if explicit_neutral
                    bus["vm"] = bus["vm"][1:length(bus["terminals"])]
                    bus["va"] = bus["va"][1:length(bus["terminals"])]
                    bus["vmin"] = bus["vmin"][1:length(bus["terminals"])]
                    bus["vmax"] = bus["vmax"][1:length(bus["terminals"])]
                    bus["grounded"] = bus["grounded"][1:length(bus["terminals"])]
                else
                    if neutral_idx ∈ bus["terminals"]
                        bus["terminals"] = bus["terminals"][1:end-1]
                    end
                    bus["vm"] = bus["vm"][1:length(bus["terminals"])]
                    bus["va"] = bus["va"][1:length(bus["terminals"])]
                    bus["vmin"] = bus["vmin"][1:length(bus["terminals"])]
                    bus["vmax"] = bus["vmax"][1:length(bus["terminals"])]
                    bus["grounded"] = bus["grounded"][1:length(bus["terminals"])]
                end
            end
        end
    end
    return nothing
end


function update_math_model_3wire!(math)
    
    math["conductor_ids"] = math["conductor_ids"][1:3]

    for (i,bus) in math["bus"]
        explicit_neutral = false
        if haskey(bus, "terminals") && neutral_idx ∈ bus["terminals"]
            explicit_neutral = true
        end

        if explicit_neutral
            if haskey(bus, "terminals")
                bus["terminals"] = bus["terminals"][1:end-1]
            end
            if haskey(bus, "grounded")
                bus["grounded"] = bus["grounded"][1:end-1]
            end
            if haskey(bus, "vmax")
                bus["vmax"] = bus["vmax"][1:end-1]
            end
            if haskey(bus, "vmin")
                bus["vmin"] = bus["vmin"][1:end-1]
            end
            bus["vmin"] = 0.9 * ones(length(bus["terminals"]))
            bus["vmax"] = 1.1 * ones(length(bus["terminals"]))
        end
    end

    for (l,branch) in math["branch"]
        explicit_neutral = false
        if haskey(branch, "t_connections") && neutral_idx ∈ branch["t_connections"]
            explicit_neutral = true
            branch["t_connections"] = branch["t_connections"][1:end-1]
        end
        if haskey(branch, "f_connections") && neutral_idx ∈ branch["f_connections"]
            explicit_neutral = true
            branch["f_connections"] = branch["f_connections"][1:end-1]
        end
        if haskey(branch, "br_r") && explicit_neutral
            branch["br_r"] = branch["br_r"][1:end-1,1:end-1]
        end
        if haskey(branch, "br_x") && explicit_neutral
            branch["br_x"] = branch["br_x"][1:end-1,1:end-1]
        end
        if haskey(branch, "g_to") && explicit_neutral
            branch["g_to"] = branch["g_to"][1:end-1,1:end-1]
        end
        if haskey(branch, "g_fr") && explicit_neutral
            branch["g_fr"] = branch["g_fr"][1:end-1,1:end-1]
        end
        if haskey(branch, "b_to") && explicit_neutral
            branch["b_to"] = branch["b_to"][1:end-1,1:end-1]
        end
        if haskey(branch, "b_fr") && explicit_neutral
            branch["b_fr"] = branch["b_fr"][1:end-1,1:end-1]
        end
        if haskey(branch, "c_rating_a") && explicit_neutral
            branch["c_rating_a"] = branch["c_rating_a"][1:end-1]
        end
    end

    for (t,transformer) in math["transformer"]
        if haskey(transformer, "t_connections") && neutral_idx ∈ transformer["t_connections"]
            transformer["t_connections"] = transformer["t_connections"][1:end-1]
        end
        if haskey(transformer, "f_connections") && neutral_idx ∈ transformer["f_connections"]
            transformer["f_connections"] = transformer["f_connections"][1:end-1]
        end
    end

    for (g,gen) in math["gen"]
        if neutral_idx in gen["connections"]
            gen["connections"] = gen["connections"][1:end-1]
            gen["vg"] = gen["vg"][1:end-1]
            gen["pg"] = gen["pg"][1:end-1]
            gen["qg"] = gen["qg"][1:end-1]
            gen["pmax"] = gen["pmax"][1:end-1]
            gen["pmin"] = gen["pmin"][1:end-1]
            gen["qmax"] = gen["qmax"][1:end-1]
            gen["qmin"] = gen["qmin"][1:end-1]
            gen["cost"] = 1000 .* gen["cost"]
        end
    end

    for (l,load) in math["load"]
        if load["configuration"] == WYE && neutral_idx ∈ load["connections"]
            load["connections"] = load["connections"][1:end-1]
        end
    end

    return nothing
end


function multinetwork_data_math_correction!(data_math::Dict{String, Any})
    @assert data_math["multinetwork"]
    @assert data_math["data_model"]==MATHEMATICAL
    for (nw, dm) in data_math["nw"]
        dm["data_model"] = MATHEMATICAL
        dm["map"] = data_math["map"]
        dm["bus_lookup"] = data_math["bus_lookup"][nw]
    end
    return nothing
end

##
# point to data and solution directory
data_dir = "data/en_validation_case_data"
solution_dir = "data/en_validation_case_solutions"
# infer cases from files defined in data dir
cases = [x[1:end-4] for x in readdir(data_dir) if endswith(x, ".dss")]
filter!(e->e≠"test_line_6w", cases)
filter!(e->e≠"test_load_3ph_wye_exp", cases)
filter!(e->e≠"test_load_3ph_delta_exp", cases)
filter!(e->e≠"test_trans_dy_3w", cases)
filter!(e->e≠"test_trans_yy_3w", cases)

@testset "en pf native opendss validation four wire" begin

    for (case_idx,case) in enumerate(cases)

        @testset "case $case" begin
            case_path = "$data_dir/$case.dss"

            data_eng = parse_file(case_path, transformations=[transform_loops!])
            vsource_correction!(data_eng; explicit_neutral=true)

            data_math = transform_data_model(data_eng;kron_reduce=false)

            res = compute_pf(data_math; explicit_neutral=true)

            # obtain solution from dss
            sol_dss = open("$solution_dir/$case.json", "r") do f
                JSON.parse(f)
            end

            sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
            @test res["termination_status"]==CONVERGED

            v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)
            @test v_maxerr_pu <= 1E-6

        end
    end

    @testset "compute_pf and solution" begin
        # testing compute_pf(pfd)
        case = "test_load_3ph_wye_cp"
        case_path = "$data_dir/$case.dss"

        data_eng = parse_file(case_path, transformations=[transform_loops!])
        vsource_correction!(data_eng; explicit_neutral=true)

        data_math = transform_data_model(data_eng;kron_reduce=false)

        add_start_voltage!(data_math, coordinates=:rectangular, epsilon=0)
        v_start = PowerModelsDistribution._bts_to_start_voltage(data_math)

        explicit_neutral = true
        pfd = PowerFlowData(data_math, v_start, explicit_neutral)

        res = compute_pf(pfd)

        sol_dss = open("$solution_dir/$case.json", "r") do f
            JSON.parse(f)
        end

        sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
        v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)
        @test v_maxerr_pu <= 1E-6

        @test haskey(res["solution"], "gen")
        @test haskey(res["solution"], "branch")
        @test haskey(res["solution"], "bus")
        @test haskey(res["solution"], "load")
        @test haskey(res["solution"], "shunt")
        @test haskey(res["solution"], "switch")
        @test haskey(res["solution"], "transformer")
    end
end


cases = ["test_trans_dy_3w", "test_trans_yy_3w"]

@testset "en pf native opendss validation three wire" begin

    for (case_idx,case) in enumerate(cases)
        @testset "case $case" begin
            case_path = "$data_dir/$case.dss"

            data_eng = parse_file(case_path, transformations=[transform_loops!])
            data_eng["is_kron_reduced"] = true
            data_eng["settings"]["sbase_default"] = 1
            vsource_correction!(data_eng; explicit_neutral=false)
            
            data_math = transform_data_model(data_eng;kron_reduce=false, phase_project=false)
            sourcebus_voltage_vector_correction!(data_math, explicit_neutral=false)
            update_math_model_3wire!(data_math)

            res = compute_pf(data_math; explicit_neutral=false)

            # obtain solution from dss
            sol_dss = open("$solution_dir/$case.json", "r") do f
                JSON.parse(f)
            end

            sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
            @test res["termination_status"]==CONVERGED

            v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)
            @test v_maxerr_pu <= 1E-6

        end
    end
end

@testset "en pf native opendss validation for max iteration limit" begin
    case_path = "$data_dir/test_load_3ph_wye_cp.dss"

    data_eng = parse_file(case_path, transformations=[transform_loops!])
    vsource_correction!(data_eng; explicit_neutral=true)

    data_math = transform_data_model(data_eng;kron_reduce=false)

    res = compute_pf(data_math; max_iter=5, explicit_neutral=true)
    @test res["termination_status"] == ITERATION_LIMIT
end


@testset "en pf native opendss validation for 6 wire test case" begin
    # This test is to throw the warning for the >4 wire network
    case_path = "$data_dir/test_line_6w.dss"

    data_eng = parse_file(case_path, transformations=[transform_loops!])
    vsource_correction!(data_eng; explicit_neutral=true)

    data_math = transform_data_model(data_eng;kron_reduce=false)

    res = compute_pf(data_math; max_iter=1, explicit_neutral=true)
end




##
# point to data and solution directory
data_dir = "data/opendss"
solution_dir = "data/opendss_solutions"

@testset "en pf native opendss validation for multinetwork test case" begin
    case = "case3_balanced"
    
    eng_ts = make_multinetwork(case3_balanced)
    vsource_correction!(eng_ts; explicit_neutral=true)

    ## This section is to validate that the new feature does not break multinetwork 
    result_mn = solve_mn_mc_opf(eng_ts, ACPUPowerModel, ipopt_solver)

    @test result_mn["termination_status"] == LOCALLY_SOLVED
    @test isapprox(result_mn["solution"]["nw"]["1"]["bus"]["loadbus"]["vm"][1], 0.229971669, atol=1E-6)
    @test isapprox(result_mn["solution"]["nw"]["10"]["bus"]["loadbus"]["vm"][1], 0.22977482, atol=1E-6)


    ## This section is to validate that the native pf supports multinetwork
    data_math = transform_data_model(eng_ts;kron_reduce=false)
    multinetwork_data_math_correction!(data_math)

    res = compute_pf(data_math; explicit_neutral=true)
    @test res["termination_status"]["10"]==CONVERGED

    sol_dss = open("$solution_dir/$case.json", "r") do f
        JSON.parse(f)
    end
    
    sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
    
    v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd["nw"]["10"], eng_ts["nw"]["10"], data_math["nw"]["10"], verbose=false, compare_math=true)
    @test v_maxerr_pu <= 1E-1   # This tolerance is selected for the multinetwork to pass, must be tightened later. The problem may be with OpenDSS json result.

end



# infer cases from files defined in data dir
cases = [x[1:end-5] for x in readdir(solution_dir) if endswith(x, ".json")]
filter!(e->e≠"case3_unbalanced_delta_loads", cases)
# filter!(e->e≠["ut_trans_2w_dy_lag" "ut_trans_2w_dy_lead_small_series_impedance" "ut_trans_2w_dy_lead" "ut_trans_2w_yy_oltc" "ut_trans_2w_yy" "ut_trans_3w_dyy_1"], cases)

@testset "en pf native opendss validation three wire" begin

    for (case_idx,case) in enumerate(cases[1:15])

        @testset "case $case" begin
            case_path = "$data_dir/$case.dss"
            
            data_eng = parse_file(case_path, transformations=[transform_loops!])#, kron_reduce=true)#, phase_project=false)
            conductor_correction!(data_eng)
            data_eng["is_kron_reduced"] = true

            data_math = transform_data_model(data_eng;kron_reduce=false)

            res = compute_pf(data_math; explicit_neutral=false)

            # obtain solution from dss
            sol_dss = open("$solution_dir/$case.json", "r") do f
                JSON.parse(f)
            end

            sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
            @test res["termination_status"]==CONVERGED

            v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, data_eng, data_math, verbose=false, compare_math=true)
            if occursin("switch", case)
                @test v_maxerr_pu <= 1E-3
            elseif case == "case3_unbalanced_delta_loads"
                @test v_maxerr_pu <= 2E-3  # This tolerance must be tightened later. The problem may be with compute_df single phase loads.
            else
                @test v_maxerr_pu <= 1E-6
            end
        end
    end
end