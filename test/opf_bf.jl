@info "running branch-flow optimal power flow (opf_bf) tests"

@testset "test distflow formulations" begin
    @testset "test linearised distflow opf_bf" begin
        @testset "5-bus lplinubf opf_bf" begin
            mp_data = PowerModels.parse_file("../test/data/matpower/case5.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf_bf(mp_data, LPLinUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 44880; atol = 1e0)
            @test isapprox(result["solution"]["bus"]["3"]["vm"].values, 0.911466*[1,1,1]; atol = 1e-3)
        end

        @testset "5-bus independent radial identical lplinubf opf_bf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_i_r_a.m")
            result = run_mc_opf_bf(mp_data, LPLinUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 54870.0; atol = 1e-1)
            @test isapprox(result["solution"]["bus"]["3"]["vm"].values, 1.02472*[1,1,1]; atol = 1e-3)
        end

        @testset "5-bus independent radial different lplinubf opf_bf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_i_r_b.m")
            result = run_mc_opf_bf(mp_data, LPLinUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 55307.7; atol = 1e-1)
            @test isapprox(result["solution"]["bus"]["3"]["vm"].values, [0.930548, 0.930543, 0.930543]; atol = 1e-3)
        end

        @testset "3-bus balanced lplinubf opf_bf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_balanced.dss")
            sol = PMD.run_mc_opf_bf(pmd, LPLinUBFPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0183456; atol=2e-3)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00923328; atol=2e-3)
        end

        @testset "3-bus unbalanced lplinubf opf_bf" begin
            pmd = PMD.parse_file("../test/data/opendss/case3_unbalanced.dss")
            sol = PMD.run_mc_opf_bf(pmd, LPLinUBFPowerModel, ipopt_solver)

            @test sol["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(sum(sol["solution"]["gen"]["1"]["pg"] * sol["solution"]["baseMVA"]), 0.0214812; atol=2e-3)
            @test isapprox(sum(sol["solution"]["gen"]["1"]["qg"] * sol["solution"]["baseMVA"]), 0.00927263; atol=2e-3)
        end
    end

    @testset "test linearised distflow opf_bf in diagonal matrix form" begin
        @testset "5-bus lpdiagubf opf_bf" begin
            mp_data = PowerModels.parse_file("../test/data/matpower/case5.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf_bf(mp_data, LPdiagUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 44880; atol = 1e0)
        end

        @testset "5-bus independent radial identical lpdiagubf opf_bf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_i_r_a.m")
            result = run_mc_opf_bf(mp_data, LPdiagUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 54870.0; atol = 1e-1)
            @test isapprox(result["solution"]["bus"]["3"]["vm"].values, 1.02472*[1,1,1]; atol = 1e-3)
        end

        @testset "5-bus independent radial different lpdiagubf opf_bf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_i_r_b.m")
            result = run_mc_opf_bf(mp_data, LPdiagUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 55307.7; atol = 1e-1)
            @test isapprox(result["solution"]["bus"]["3"]["vm"].values, [0.930014, 0.930014, 0.930014]; atol = 1e-3)
        end
    end

    @testset "test linearised distflow opf_bf in full matrix form" begin
        @testset "5-bus lpfullubf opf_bf" begin
            mp_data = PowerModels.parse_file("../test/data/matpower/case5.m")
            PMD.make_multiconductor!(mp_data, 3)
            result = run_mc_opf_bf(mp_data, LPfullUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 44880; atol = 1e0)
        end

        @testset "5-bus independent radial identical lpfullubf opf_bf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_i_r_a.m")
            result = run_mc_opf_bf(mp_data, LPfullUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 54870.0; atol = 1e-1)
        end

        @testset "5-bus independent radial different lpfullubf opf_bf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_i_r_b.m")
            result = run_mc_opf_bf(mp_data, LPfullUBFPowerModel, ipopt_solver)

            @test result["termination_status"] == PMs.LOCALLY_SOLVED
            @test isapprox(result["objective"], 55307.7; atol = 1e-1)
            @test isapprox(result["solution"]["bus"]["3"]["vm"].values, [0.930014, 0.930014, 0.930014]; atol = 1e-3)

        end
    end

    @testset "test sdp distflow opf_bf" begin
        @testset "5-bus independent radial identical sdpubf opf_bf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_i_r_a.m")
            result = run_mc_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

            @test result["termination_status"] == PMs.OPTIMAL
            @test isapprox(result["objective"], 55451.2; atol = 2e0)
        end

        @testset "5-bus independent radial different sdpubf opf_bf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_i_r_b.m")
            result = run_mc_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

            @test result["termination_status"] == PMs.OPTIMAL
            @test isapprox(result["objective"], 56091.7; atol = 2e0)

            @test isapprox(result["solution"]["bus"]["1"]["vm"][1],   1.076184696745133; atol = 1E-3)
            @test isapprox(result["solution"]["bus"]["2"]["vm"][1],   1.063950148862639; atol = 1E-3)
            @test isapprox(result["solution"]["bus"]["3"]["vm"][1],   1.073074502462763; atol = 1E-3)
            @test isapprox(result["solution"]["bus"]["4"]["vm"][1],   1.053610154771143; atol = 1E-3)
        end

        @testset "5-bus independent meshed different sdpubf opf_bf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_i_m_b.m")
            result = run_mc_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

            @test result["termination_status"] == PMs.OPTIMAL
            @test isapprox(result["objective"], 45320.5; atol = 2e0)
        end

        @testset "5-bus coupled meshed infeasible sdpubf opf_bf" begin
            @testset "ac case" begin
                mp_data = PMD.parse_file("../test/data/matlab/case5_c_m_b.m")
                result = run_mc_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

                @test result["termination_status"] == PMs.INFEASIBLE
            end
        end

        @testset "5-bus coupled radial no shunt sdpubf opf_bf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_c_r_a.m")
            result = run_mc_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

            @test result["termination_status"] == PMs.OPTIMAL
            @test isapprox(result["objective"], 55434.8; atol = 2e1)

            for c in 1:mp_data["conductors"]
                @test isapprox(result["solution"]["gen"]["1"]["pg"][c], 0.4; atol = 1e-3)
                @test isapprox(result["solution"]["bus"]["2"]["vm"][c], 1.08620; atol = 1e-3)
            end
        end

        @testset "5-bus coupled radial shunt sdpubf opf_bf" begin
            mp_data = PMD.parse_file("../test/data/matlab/case5_c_r_b.m")
            result = run_mc_opf_bf(mp_data, SDPUBFPowerModel, scs_solver)

            @test result["termination_status"] == PMs.OPTIMAL
            @test isapprox(result["objective"], 56075.9; atol = 2e0)

            for c in 1:mp_data["conductors"]
                @test isapprox(result["solution"]["gen"]["1"]["pg"][c],  0.4; atol = 1e-3)
            end
        end
    end
end
