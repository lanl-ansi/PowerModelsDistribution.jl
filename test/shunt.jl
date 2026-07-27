@info "running matrix shunt tests"

@testset "matrix shunts ACP/ACR/IVR" begin
    data = transform_data_model(case2_mxshunt_2)

    shunt = data["shunt"]["1"]

    @test(isa(shunt["gs"], Matrix))
    @test(isa(shunt["bs"], Matrix))

    @test(all([shunt["bs"][c,d]!=0 for c in 1:3, d in 1:3 if c!=d]))

    data_diag_shunt = deepcopy(data)

    data_diag_shunt["shunt"]["1"]["bs"] = shunt["bs"].*[c==d for c in 1:3, d in 1:3]

    result_acp_diag = solve_mc_pf(data_diag_shunt, ACPUPowerModel, ipopt_solver)
    result_acp = solve_mc_pf(data, ACPUPowerModel, ipopt_solver)
    result_acr = solve_mc_pf(data, ACRUPowerModel, ipopt_solver)
    result_iv = solve_mc_pf(data, IVRUPowerModel, ipopt_solver)

    # check the results are different with only diagonal elements
    @test(!isapprox(result_acp["solution"]["bus"]["2"]["vm"], result_acp_diag["solution"]["bus"]["2"]["vm"]))
    # now check that the other non-linear formulations give the same result
    @test(isapprox(calc_vm_acr(result_acr, data, "loadbus"), result_acp["solution"]["bus"]["2"]["vm"], atol=1E-6))
    @test(isapprox(calc_vm_acr(result_iv, data, "loadbus"), result_acp["solution"]["bus"]["2"]["vm"], atol=1E-6))
end
