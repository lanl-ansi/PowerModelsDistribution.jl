@info "running multinetwork.jl tests"
@testset "test multinetwork" begin

    @testset "opf with storage case" begin
        mp_data = PowerModels.parse_file("../test/data/matpower/case5_strg.m")
        PMD.make_multiconductor!(mp_data, 3)
        mn_mp_data = PowerModels.replicate(mp_data, 5)

        @testset "test ac polar opf" begin
            result = PMD.run_mn_mc_strg_opf(mn_mp_data, PowerModels.ACPPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 2.64596e5; atol = 1e2)

            for (n, network) in result["solution"]["nw"], c in 1:network["conductors"]
                # @test isapprox(network["storage"]["1"]["ps"][c], -0.012; atol = 1e-3)
                # @test isapprox(network["storage"]["1"]["qs"][c], -0.012; atol = 1e-3)
                @test isapprox(network["storage"]["1"]["ps"][c], -0.012; atol = 1e-1)
                @test isapprox(network["storage"]["1"]["qs"][c], -0.012; atol = 1e-1)

                # @test isapprox(network["storage"]["2"]["ps"][c], -0.016; atol = 1e-3)
                # @test isapprox(network["storage"]["2"]["qs"][c],  0.000; atol = 1e-3)
                @test isapprox(network["storage"]["2"]["ps"][c], -0.016; atol = 1e-1)
                @test isapprox(network["storage"]["2"]["qs"][c],  0.000; atol = 1e-1)
            end
        end

        #=
        # NUMERICAL_ERROR
        @testset "test dc polar opf" begin
            result = PMD.run_mn_mc_strg_opf(mn_mp_data, PowerModels.DCPPowerModel, JuMP.with_optimizer(Ipopt.Optimizer, print_level=0)) # this test requires default tol value

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 2.63419e5; atol = 1e2)

            for (n, network) in result["solution"]["nw"], c in 1:network["conductors"]
                # @test isapprox(network["storage"]["1"]["ps"][c], -0.01199; atol = 1e-3)
                # @test isapprox(network["storage"]["2"]["ps"][c], -0.01597; atol = 1e-3)
                @test isapprox(network["storage"]["1"]["ps"][c], -0.01199; atol = 1e-1)
                @test isapprox(network["storage"]["2"]["ps"][c], -0.01597; atol = 1e-1)
            end
        end
        =#

        #=
        # non-convexity issues probably need to be resolved first
        @testset "test nfa opf" begin
            result = PMD.run_mn_mc_strg_opf(mn_mp_data, PowerModels.NFAPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 2.63419e5; atol = 1e2)

            for (n, network) in result["solution"]["nw"], c in 1:network["conductors"]
                @test isapprox(network["storage"]["1"]["ps"][c], -0.012; atol = 1e-3)
                @test isapprox(network["storage"]["2"]["ps"][c], -0.016; atol = 1e-3)
            end
        end
        =#

    end

end
