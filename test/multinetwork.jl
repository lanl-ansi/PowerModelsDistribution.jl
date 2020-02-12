@info "running multinetwork tests"

@testset "test multinetwork" begin
    # TODO: Add new Multinetwork tests after support for LoadShapes
#     @testset "5-bus storage multinetwork acp opf_strg" begin
#         mp_data = PowerModels.parse_file("../test/data/matpower/case5_strg.m")
#         PMD.make_multiconductor!(mp_data, 3)
#         mn_mp_data = PowerModels.replicate(mp_data, 5)

#         result = PMD._run_mn_mc_opf_strg(mn_mp_data, PowerModels.ACPPowerModel, ipopt_solver)

#         @test result["termination_status"] == PMs.LOCALLY_SOLVED
#         @test isapprox(result["objective"], 2.64596e5; atol = 1e2)

#         for (n, network) in result["solution"]["nw"], c in 1:network["conductors"]
#             # @test isapprox(network["storage"]["1"]["ps"][c], -0.012; atol = 1e-3)
#             # @test isapprox(network["storage"]["1"]["qs"][c], -0.012; atol = 1e-3)
#             @test isapprox(network["storage"]["1"]["ps"][c], -0.012; atol = 1e-1)
#             # @test isapprox(network["storage"]["1"]["qs"][c], -0.012; atol = 1e-1)

#             # @test isapprox(network["storage"]["2"]["ps"][c], -0.016; atol = 1e-3)
#             # @test isapprox(network["storage"]["2"]["qs"][c],  0.000; atol = 1e-3)
#             @test isapprox(network["storage"]["2"]["ps"][c], -0.016; atol = 1e-1)
#             @test isapprox(network["storage"]["2"]["qs"][c],  0.000; atol = 1e-1)
#         end
#     end
end
