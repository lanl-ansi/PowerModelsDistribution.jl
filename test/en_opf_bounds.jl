@info "running explicit neutral opf bound tests"

function calc_sol_pmd(data_math, form; optimizer=ipopt_solver)
    pm  = instantiate_mc_model(data_math, form, build_mc_opf)
    res = optimize_model!(pm, optimizer=optimizer)
    @test res["termination_status"] == LOCALLY_SOLVED
    sol_pmd = transform_solution(res["solution"], data_math, make_si=true)
    return sol_pmd
end



@testset "en opf bounds" begin
    @testset "branch current magnitude bound" begin
        cm_ub = [6:-1:3...]
        data_eng = deepcopy(test_gen_3ph_wye)
        data_eng["settings"]["sbase_default"] = 1.0
        data_eng["line"]["line1"]["cm_ub"] = cm_ub
        data_math = transform_data_model(data_eng, multinetwork=false, kron_reduced=false, phase_projected=false)
        add_start_vrvi!(data_math)
        gen_pv = data_math["gen"]["1"]
        @assert gen_pv["name"]=="pv1"
        gen_pv["cost"] *= 0.5
        gen_pv["pmax"] = fill(Inf, 3)
        # scale up objective to prevent feasibility issues
        for (_,gen) in data_math["gen"]
            gen["cost"] *= 1E3
        end

        # IVRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRENPowerModel)
        c_to = sol_pmd["line"]["line1"]["cr_to"]+im*sol_pmd["line"]["line1"]["ci_to"]
        @test all(isapprox.(abs.(c_to[1:3]), cm_ub[1:3], rtol=0.005))

        # IVRQuadraticENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRQuadraticENPowerModel)
        c_to = sol_pmd["line"]["line1"]["cr_to"]+im*sol_pmd["line"]["line1"]["ci_to"]
        @test all(isapprox.(abs.(c_to[1:3]), cm_ub[1:3], rtol=0.005))

        # IVRReducedQuadraticENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRReducedQuadraticENPowerModel)
        c_to = sol_pmd["line"]["line1"]["cr_to"]+im*sol_pmd["line"]["line1"]["ci_to"]
        @test all(isapprox.(abs.(c_to[1:3]), cm_ub[1:3], rtol=0.005))

        # ACRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, ACRENPowerModel, optimizer=ipopt_solver)
        s_to = sol_pmd["line"]["line1"]["pt"]+im*sol_pmd["line"]["line1"]["qt"]
        v_to = sol_pmd["bus"]["b2"]["vr"]+im*sol_pmd["bus"]["b2"]["vi"]
        c_to = conj.(s_to./v_to)
        @test all(isapprox.(abs.(c_to[1:3]), cm_ub[1:3], rtol=0.005))
    end

    @testset "branch power magnitude bound" begin
        sm_ub = [90:-20:30...]
        data_eng = deepcopy(test_gen_3ph_wye)
        data_eng["settings"]["sbase_default"] = 1.0
        data_eng["line"]["line1"]["sm_ub"] = sm_ub
        data_math = transform_data_model(data_eng, multinetwork=false, kron_reduced=false, phase_projected=false)
        add_start_vrvi!(data_math)
        gen_pv = data_math["gen"]["1"]
        @assert gen_pv["name"]=="pv1"
        gen_pv["cost"] *= 0.5
        gen_pv["pmax"] = fill(Inf, 3)
        # scale up objective to prevent feasibility issues
        for (_,gen) in data_math["gen"]
            gen["cost"] *= 1E3
        end

        # IVRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRENPowerModel)
        s_to = sol_pmd["line"]["line1"]["pt"]+im*sol_pmd["line"]["line1"]["qt"]
        @test all(isapprox.(abs.(s_to[1:3]), sm_ub[1:3], rtol=0.005))

        # IVRQuadraticENPowerModel does not implement sm_ub

        # IVRReducedQuadraticENPowerModel does not implement sm_ub

        # ACRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, ACRENPowerModel)
        s_to = sol_pmd["line"]["line1"]["pt"]+im*sol_pmd["line"]["line1"]["qt"]
        @test all(isapprox.(abs.(s_to[1:3]), sm_ub[1:3], rtol=0.005))
    end

    @testset "switch current magnitude bound" begin
        cm_ub = [10:-2:4...]
        data_eng = deepcopy(test_switch)
        add_bus_absolute_vbounds!(data_eng) # ACR needs this to prevent voltage collapse
        # copy in solar from test_gen_3ph_wye.dss
        data_eng["solar"] = deepcopy(test_gen_3ph_wye["solar"])
        data_eng["settings"]["sbase_default"] = 500.0
        data_eng["switch"]["switch"]["cm_ub"] = cm_ub
        data_math = transform_data_model(data_eng, multinetwork=false, kron_reduced=false, phase_projected=false)
        add_start_vrvi!(data_math)
        gen_pv = data_math["gen"]["1"]
        @assert gen_pv["name"]=="pv1"
        gen_pv["cost"] *= 0.5
        gen_pv["pmax"] = fill(Inf, 3)
        # scale up objective to prevent feasibility issues
        for (_,gen) in data_math["gen"]
            gen["cost"] *= 1E4
        end

        # IVRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRENPowerModel)
        c_to = sol_pmd["switch"]["switch"]["cr_to"]+im*sol_pmd["switch"]["switch"]["ci_to"]
        @test all(isapprox.(abs.(c_to[1:3]), cm_ub[1:3], rtol=0.005))

        # IVRQuadraticENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRQuadraticENPowerModel)
        c_to = sol_pmd["switch"]["switch"]["cr_to"]+im*sol_pmd["switch"]["switch"]["ci_to"]
        @test all(isapprox.(abs.(c_to[1:3]), cm_ub[1:3], rtol=0.005))

        # IVRReducedQuadraticENPowerModel shares switch implementation with IVRQuadraticENPowerModel,
        # so no explicit test needed

        # ACRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, ACRENPowerModel)
        s_to = sol_pmd["switch"]["switch"]["pt"]+im*sol_pmd["switch"]["switch"]["qt"]
        v_to = sol_pmd["bus"]["x2"]["vr"]+im*sol_pmd["bus"]["x2"]["vi"]
        c_to = conj.(s_to./v_to)
        @test all(isapprox.(abs.(c_to[1:3]), cm_ub[1:3], rtol=0.01))
    end

    @testset "branch power magnitude bound" begin
        sm_ub = [50:-10:20...]
        data_eng = deepcopy(test_switch)
        add_bus_absolute_vbounds!(data_eng) # ACR needs this to prevent voltage collapse
        # copy in solar from test_gen_3ph_wye.dss
        data_eng["solar"] = deepcopy(test_gen_3ph_wye["solar"])
        data_eng["settings"]["sbase_default"] = 500.0
        data_eng["switch"]["switch"]["sm_ub"] = sm_ub
        data_math = transform_data_model(data_eng, multinetwork=false, kron_reduced=false, phase_projected=false)
        add_start_vrvi!(data_math)
        gen_pv = data_math["gen"]["1"]
        @assert gen_pv["name"]=="pv1"
        gen_pv["cost"] *= 0.5
        gen_pv["pmax"] = fill(Inf, 3)
        # scale up objective to prevent feasibility issues
        for (_,gen) in data_math["gen"]
            gen["cost"] *= 1E3
        end

        # IVRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRENPowerModel, optimizer=optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
        c_to = sol_pmd["switch"]["switch"]["cr_to"]+im*sol_pmd["switch"]["switch"]["ci_to"]
        v_to = sol_pmd["bus"]["x2"]["vr"]+im*sol_pmd["bus"]["x2"]["vi"]
        s_to = v_to.*conj.(c_to)
        @test all(isapprox.(abs.(s_to[1:3]), sm_ub[1:3], rtol=0.01))

        # IVRQuadraticENPowerModel does not implement sm_ub

        # IVRReducedQuadraticENPowerModel does not implement sm_ub

        # ACRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, ACRENPowerModel)
        s_to = sol_pmd["switch"]["switch"]["pt"]+im*sol_pmd["switch"]["switch"]["qt"]
        v_to = sol_pmd["bus"]["x2"]["vr"]+im*sol_pmd["bus"]["x2"]["vi"]
        @test all(isapprox.(abs.(s_to[1:3]), sm_ub[1:3], rtol=0.001))
    end

    @testset "bus absolute voltage magnitude upper bound" begin
        vm_ub = [0.300, 0.310, 0.320, 0.010]
        data_eng = deepcopy(test_gen_3ph_wye)
        data_eng["settings"]["sbase_default"] = 500.0
        data_eng["bus"]["b2"]["vm_ub"] = vm_ub
        data_math = transform_data_model(data_eng, multinetwork=false, kron_reduced=false, phase_projected=false)
        add_start_vrvi!(data_math)
        gen_pv = data_math["gen"]["1"]
        @assert gen_pv["name"]=="pv1"
        gen_pv["cost"] *= 0.5
        gen_pv["pmax"] = fill(Inf, 3)
        # scale up objective to prevent feasibility issues
        for (_,gen) in data_math["gen"]
            gen["cost"] *= 1E3
        end

        # IVRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRENPowerModel)
        v_b2 = sol_pmd["bus"]["b2"]["vr"]+im*sol_pmd["bus"]["b2"]["vi"]
        @test all(isapprox.(abs.(v_b2[1:3]), vm_ub[1:3], rtol=0.005))

        # IVRQuadraticENPowerModel, IVRReducedQuadraticENPowerModel and ACRENPowerModel
        # share switch implementation with IVRQuadraticENPowerModel,
        # so no explicit test needed
    end

    @testset "bus absolute voltage magnitude lower bound" begin
        vm_lb = [0.18, 0.190, 0.170, 0.0]
        data_eng = deepcopy(test_gen_3ph_wye)
        data_eng["settings"]["sbase_default"] = 500.0
        data_eng["bus"]["b2"]["vm_lb"] = vm_lb
        data_math = transform_data_model(data_eng, multinetwork=false, kron_reduced=false, phase_projected=false)
        add_start_vrvi!(data_math)
        gen_pv = data_math["gen"]["1"]
        @assert gen_pv["name"]=="pv1"
        gen_pv["cost"] *= 2.0 # reverse power flow
        gen_pv["pmax"] = fill( Inf, 3)
        gen_pv["pmin"] = fill(-Inf, 3)
        gen_pv["qmin"] = gen_pv["qmax"] = fill(0.0, 3)
        # scale up objective to prevent feasibility issues
        for (_,gen) in data_math["gen"]
            gen["cost"] *= 1E3
        end

        # IVRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRENPowerModel)
        v_b2 = sol_pmd["bus"]["b2"]["vr"]+im*sol_pmd["bus"]["b2"]["vi"]
        @test any(isapprox.(abs.(v_b2[1:3]), vm_lb[1:3], rtol=0.005))

        # IVRQuadraticENPowerModel, IVRReducedQuadraticENPowerModel and ACRENPowerModel
        # share switch implementation with IVRQuadraticENPowerModel,
        # so no explicit test needed
    end

    @testset "bus pairwise voltage magnitude upper bound" begin
        vm_pair_ub = [(1,4,0.250)]
        data_eng = deepcopy(test_gen_1ph_wye)
        data_eng["settings"]["sbase_default"] = 500.0
        data_eng["bus"]["b2"]["vm_pair_ub"] = vm_pair_ub
        data_math = transform_data_model(data_eng, multinetwork=false, kron_reduced=false, phase_projected=false)
        add_start_vrvi!(data_math)
        gen_pv = data_math["gen"]["1"]
        @assert gen_pv["name"]=="pv1"
        gen_pv["cost"] *= 0.5
        gen_pv["pmax"] = fill(Inf, 3)
        # scale up objective to prevent feasibility issues
        for (_,gen) in data_math["gen"]
            gen["cost"] *= 1E3
        end

        # IVRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRENPowerModel)
        v_b2 = sol_pmd["bus"]["b2"]["vr"]+im*sol_pmd["bus"]["b2"]["vi"]
        c, d, lb = vm_pair_ub[1]
        @test all(isapprox.(abs(v_b2[c]-v_b2[d]), lb, rtol=0.005))

        # IVRQuadraticENPowerModel, IVRReducedQuadraticENPowerModel and ACRENPowerModel
        # share switch implementation with IVRQuadraticENPowerModel,
        # so no explicit test needed
    end

    @testset "bus pairwise voltage magnitude lower bound" begin
        vm_pair_lb = [(1,4,0.150)]
        data_eng = deepcopy(test_gen_1ph_wye)
        data_eng["settings"]["sbase_default"] = 500.0
        data_eng["bus"]["b2"]["vm_pair_lb"] = vm_pair_lb
        data_math = transform_data_model(data_eng, multinetwork=false, kron_reduced=false, phase_projected=false)
        add_start_vrvi!(data_math)
        gen_pv = data_math["gen"]["1"]
        @assert gen_pv["name"]=="pv1"
        gen_pv["cost"] *= 10.0 # reverse power flow
        gen_pv["pmax"] = [Inf]
        gen_pv["qmax"] = [0.0]
        # scale up objective to prevent feasibility issues
        for (_,gen) in data_math["gen"]
            gen["cost"] *= 1E3
        end

        # IVRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRENPowerModel)
        v_b2 = sol_pmd["bus"]["b2"]["vr"]+im*sol_pmd["bus"]["b2"]["vi"]
        c, d, lb = vm_pair_lb[1]
        @test all(isapprox.(abs(v_b2[c]-v_b2[d]), lb, rtol=0.005))

        # IVRQuadraticENPowerModel, IVRReducedQuadraticENPowerModel and ACRENPowerModel
        # share switch implementation with IVRQuadraticENPowerModel,
        # so no explicit test needed
    end

    @testset "transformer total power magnitude upper bound" begin
        sm_ub = 100.0
        data_eng = deepcopy(test_trans_dy)
        # copy in solar from test_gen_3ph_wye.dss
        data_eng["solar"] = deepcopy(test_gen_3ph_wye["solar"])
        data_eng["settings"]["sbase_default"] = 500.0
        data_eng["transformer"]["transformer1"]["sm_ub"] = sm_ub
        data_math = transform_data_model(data_eng, multinetwork=false, kron_reduced=false, phase_projected=false)
        add_start_vrvi!(data_math)
        gen_pv = data_math["gen"]["1"]
        @assert gen_pv["name"]=="pv1"
        gen_pv["cost"] *= 0.5
        gen_pv["pmax"] = fill(Inf, 3)
        # scale up objective to prevent feasibility issues
        for (_,gen) in data_math["gen"]
            gen["cost"] *= 1E3
        end

        # IVRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRENPowerModel)
        s_to = sol_pmd["transformer"]["transformer1"]["p"][2]+im*sol_pmd["transformer"]["transformer1"]["q"][2]
        @assert isapprox(abs(sum(s_to)), sm_ub, rtol=0.005)

        # IVRQuadraticENPowerModel
        sol_pmd = calc_sol_pmd(data_math, IVRQuadraticENPowerModel)
        s_to = sol_pmd["transformer"]["transformer1"]["p"][2]+im*sol_pmd["transformer"]["transformer1"]["q"][2]
        @assert isapprox(abs(sum(s_to)), sm_ub, rtol=0.005)

        # IVRReducedQuadraticENPowerModel shares switch implementation with IVRQuadraticENPowerModel,
        # so no explicit test needed

        # ACRENPowerModel
        sol_pmd = calc_sol_pmd(data_math, ACRENPowerModel)
        s_to = sol_pmd["transformer"]["transformer1"]["p"][2]+im*sol_pmd["transformer"]["transformer1"]["q"][2]
        @assert isapprox(abs(sum(s_to)), sm_ub, rtol=0.005)
    end
end