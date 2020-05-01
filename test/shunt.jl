@info "running matrix shunt tests"

@testset "matrix shunts ACP/ACR/IVR" begin
    data = parse_file("data/opendss/case_mxshunt_2.dss"; data_model=MATHEMATICAL)
    shunt = data["shunt"]["1"]
    @test(isa(shunt["gs"], Matrix))
    @test(isa(shunt["bs"], Matrix))
    @test(all([shunt["bs"][c,d]!=0 for c in 1:3, d in 1:3 if c!=d]))

    data_diag_shunt = deepcopy(data)
    data_diag_shunt["shunt"]["1"]["bs"] = shunt["bs"].*[c==d for c in 1:3, d in 1:3]

    sol_acp_diag = run_mc_pf(data_diag_shunt, ACPPowerModel, ipopt_solver)
    sol_acp = run_mc_pf(data, ACPPowerModel, ipopt_solver)
    sol_acr = run_mc_pf(data, ACRPowerModel, ipopt_solver)
    sol_iv = run_mc_pf(data, IVRPowerModel, ipopt_solver)

    # check the results are different with only diagonal elements
    @test(!isapprox(sol_acp["solution"]["bus"]["2"]["vm"], sol_acp_diag["solution"]["bus"]["2"]["vm"]))
    # now check that the other non-linear formulations give the same result
    @test(isapprox(calc_vm_acr(sol_acr, data, "loadbus"), sol_acp["solution"]["bus"]["2"]["vm"], atol=1E-6))
    @test(isapprox(calc_vm_acr(sol_iv, data, "loadbus"), sol_acp["solution"]["bus"]["2"]["vm"], atol=1E-6))
end
