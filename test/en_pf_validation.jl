@info "running explicit neutral power flow tests"

"The solar parsing is a bit off, this method corrects for that."
function pv1_correction!(data_eng)
    pv1 = data_eng["solar"]["pv1"]
    pv1["pg_lb"] = pv1["pg_ub"] = pv1["pg"]
    pv1["qg_lb"] = pv1["qg_ub"] = pv1["qg"]
end

"""
The ACR formulation returns non-physical solutions when the neutral voltage is not bounded below.
This method uses the ODD solution to find valid lower bounds, so the formulation can be validated.
"""
function add_neutral_lb_from_soldss(data_eng, sol_dss)
    de = deepcopy(data_eng)
    for (id, sol_bus) in sol_dss["bus"]
        eng_bus = de["bus"][id]
        ts = eng_bus["terminals"]
        for t in [t for t in ts if t>3 && t in keys(sol_bus["vm"])]
            idx = findfirst(ts.==t)
            eng_bus["vm_lb"] = haskey(eng_bus, "vm_lb") ? eng_bus["vm_lb"] : fill(0.0, length(ts))
            eng_bus["vm_lb"][idx] = max(eng_bus["vm_lb"][idx], sol_bus["vm"][t]*1E-3*0.9)
        end
    end
    return de
end

"Compares a PMD and OpenDSS solution, and returns the largest difference in voltage profile in per unit."
function compare_sol_dss_pmd(sol_dss::Dict{String,Any}, sol_pmd::Dict{String,Any}, data_eng::Dict{String,Any}, data_math::Dict{String,Any}; compare_math=false, verbose=true, floating_buses=[], skip_buses=[], v_err_print_tol=1E-6)
    max_v_err_pu = 0.0

    # voltage base for ENGINEERING buses in [V]
    vbase = Dict(id=>data_math["bus"]["$ind"]["vbase"]*data_math["settings"]["voltage_scale_factor"] for (id,ind) in data_math["bus_lookup"])

    buses_intersected = intersect(keys(sol_dss["bus"]), keys(sol_pmd["bus"]))
    for id in setdiff(buses_intersected, skip_buses)
        pmd_bus = sol_pmd["bus"][id]
        dss_bus = sol_dss["bus"][id]

        terminals = data_eng["bus"][id]["terminals"]
        if compare_math
            ts = [t for t in string.(terminals) if haskey(dss_bus["vm"], t)]
            v_dss = [dss_bus["vm"][t]*exp(im*dss_bus["va"][t]) for t in ts]
            # convert to V instead of usual kV
            v_pmd = [pmd_bus["vm"][t]*exp(im*deg2rad(pmd_bus["va"][t]))*data_eng["settings"]["voltage_scale_factor"] for t in ts]
        else
            ts = [t for t in terminals if haskey(dss_bus["vm"], t)]
            v_dss = [dss_bus["vm"][t]*exp(im*dss_bus["va"][t]) for t in ts]
            # convert to V instead of usual kV
            v_pmd = [(pmd_bus["vr"][idx]+im*pmd_bus["vi"][idx])*data_eng["settings"]["voltage_scale_factor"] for (idx,t) in enumerate(ts)]
        end
        
        # convert to pu
        v_dss_pu = v_dss/vbase[id]
        v_pmd_pu = v_pmd/vbase[id]

        # convert to diffs if floating
        N = length(v_dss)
        if id in floating_buses && N>1
            v_dss_pu = [v_dss_pu[i]-v_dss_pu[j] for i in 1:N for j in i:N if i!=j]/vbase[id]
            v_pmd_pu = [v_pmd_pu[i]-v_pmd_pu[j] for i in 1:N for j in i:N if i!=j]/vbase[id]
            labels = ["$(ts[i])-$(ts[j])" for i in 1:N for j in i:N if i!=j]
        else
            labels = string.(ts)
        end

        for i in 1:length(v_pmd_pu)
            v_err_pu = abs.(v_dss_pu[i]-v_pmd_pu[i]); max_v_err_pu = max(max_v_err_pu, v_err_pu)

            if v_err_pu>v_err_print_tol && verbose
                println("terminal $id.$(labels[i])")
                println("\t |U| dss: $(abs(v_dss_pu[i]))")
                println("\t     pmd: $(abs(v_pmd_pu[i]))")
                println("\t  ∠U dss:  $(angle(v_dss_pu[i]))")
                println("\t     pmd:  $(angle(v_pmd_pu[i]))")
            end
        end
    end

    return max_v_err_pu
end

# formulations to check
forms = [IVRENPowerModel, IVRReducedENPowerModel, IVRQuadraticENPowerModel, IVRReducedQuadraticENPowerModel, ACRENPowerModel]
# point to data and solution directory
data_dir = "data/en_validation_case_data"
solution_dir = "data/en_validation_case_solutions"
# infer cases from files defined in data dir
cases = [x[1:end-4] for x in readdir(data_dir) if endswith(x, ".dss")]
# list required transformations per case
case_transformations = Dict(
    "test_gen_3ph_wye" => [remove_all_bounds!, pv1_correction!],
    "test_gen_3ph_delta" => [remove_all_bounds!, pv1_correction!],
    "test_gen_1ph_wye" => [remove_all_bounds!, pv1_correction!],
    "test_gen_1ph_wye_debug" => [remove_all_bounds!, pv1_correction!],
    "test_gen_1ph_delta" => [remove_all_bounds!, pv1_correction!],
    "test_line_6w" => [remove_all_bounds!, transform_loops!],
    "test_grounding" => [remove_all_bounds!, transform_loops!],
    "test_trans_dy" => [
        remove_all_bounds!, transform_loops!,
        x->add_bus_absolute_vbounds!(x, phase_lb_pu=0.85), # ACR needs some additional help to converge
    ],
)

@testset "en pf opendss validation" begin

    for (case_idx,case) in enumerate(cases)

        @testset "case $case" begin
            case_path = "$data_dir/$case.dss"

            transformations = haskey(case_transformations, case) ? case_transformations[case] : [remove_all_bounds!]
            data_eng = parse_file(case_path, transformations=transformations)

            # obtain solution from dss
            sol_dss = open("$solution_dir/$case.json", "r") do f
                JSON.parse(f)
            end

            # add lb on neutrals to prevent issues with ACR formulations
            data_eng_lb = add_neutral_lb_from_soldss(data_eng, sol_dss)

            # apply data model transformation and add voltage initialization
            data_math    = transform_data_model(data_eng, multinetwork=false, kron_reduce=false, phase_project=false)
            add_start_vrvi!(data_math)
            data_math_lb = transform_data_model(data_eng_lb, multinetwork=false, kron_reduce=false, phase_project=false)
            add_start_vrvi!(data_math_lb)

            for form in forms
                # for ACR formulations, use version with neutral bounds
                dm = form <: AbstractExplicitNeutralACRModel ? data_math_lb : data_math
                de = form <: AbstractExplicitNeutralACRModel ? data_eng_lb : data_eng

                pm  = instantiate_mc_model(dm, form, build_mc_opf)
                res = optimize_model!(pm, optimizer=ipopt_solver)
                sol_pmd = transform_solution(res["solution"], dm, make_si=true)
                # @assert res["termination_status"]==LOCALLY_SOLVED

                v_maxerr_pu = compare_sol_dss_pmd(sol_dss, sol_pmd, de, dm, verbose=false)
                @test v_maxerr_pu <= 1E-6
            end
        end
    end
end

