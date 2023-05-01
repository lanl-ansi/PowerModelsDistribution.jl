@info "running explicit neutral power flow tests"

# formulations to check
forms = [IVRENPowerModel, IVRReducedENPowerModel, IVRQuadraticENPowerModel, IVRReducedQuadraticENPowerModel, ACRENPowerModel]
# point to data and solution directory
data_dir = "data/en_validation_case_data"
solution_dir = "data/en_validation_case_solutions"
# infer cases from files defined in data dir
cases = [x[1:end-4] for x in readdir(data_dir) if endswith(x, ".dss")]
filter!(e->e≠"test_trans_yy", cases)
filter!(e->e≠"test_trans_dy_3w", cases)
filter!(e->e≠"test_trans_yy_3w", cases)
filter!(e -> e ≠ "case3_balanced_battery_1ph", cases)
filter!(e -> e ≠ "case3_balanced_battery_3ph", cases)
filter!(e -> e ≠ "case3_balanced_battery_3ph_EN", cases)
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

