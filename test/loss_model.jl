@testset "_build_loss_model! transformer winding placement" begin

    # ------------------------------------------------------------------
    # Starting point:
    #   data_eng should already be the parsed engineering model for TX1.
    #
    # Example:
    #   data_eng = PMD.parse_opendss("path/to/ut_trans.dss")
    #
    # The transformer should have:
    #   vm_nom = [11.0, 4.0, 0.4]
    #   sm_nom = [500.0, 500.0, 500.0]
    #   xsc    = [0.04, 0.06, 0.08]
    # ------------------------------------------------------------------

    file = "../test/data/opendss/ut_trans_3w_dyy_basetest.dss"
    data_eng = open(file) do io
       parse_opendss(io)
       end
       
    data_eng["transformer"]["tx1"]["xsc"] = [0.04, 0.06, 0.08]
    tx0 = data_eng["transformer"]["tx1"]

    @test tx0["vm_nom"] ≈ [11.0, 4.0, 0.4]
    @test tx0["sm_nom"] ≈ [500.0, 500.0, 500.0]
    @test collect(skipmissing(tx0["xsc"])) ≈ [0.04, 0.06, 0.08]

    # Get all mathematical branches derived from TX1.
    function tx1_branches(data_math)
        return [
            (id, br)
            for (id, br) in data_math["branch"]
            if occursin(
                "transformer.tx1",
                lowercase(string(get(br, "source_id", "")))
            )
        ]
    end

    # Get all mathematical ideal transformers derived from TX1.
    function tx1_transformers(data_math)
        return [
            (id, tr)
            for (id, tr) in data_math["transformer"]
            if occursin(
                "transformer.tx1",
                lowercase(string(get(tr, "source_id", "")))
            )
        ]
    end

    # Maximum absolute series resistance/reactance in a branch matrix.
    rmag(br) = maximum(abs, br["br_r"])
    xmag(br) = maximum(abs, br["br_x"])

    # Cases use ENGINEERING-MODEL per-unit winding resistances.
    #
    # With equal 500 kVA ratings, the values reaching _build_loss_model!
    # should be:
    #
    #   rw = 0.01 -> r_s = 0.02
    #   rw = 0.02 -> r_s = 0.04
    #   rw = 0.03 -> r_s = 0.06
    #
    cases = [
        (
            name = "winding 1 resistance only",
            rw = [0.01, 0.0, 0.0],
            expected_r = 0.02,
            winding = 1,
        ),
        (
            name = "winding 2 resistance only",
            rw = [0.0, 0.02, 0.0],
            expected_r = 0.04,
            winding = 2,
        ),
        (
            name = "winding 3 resistance only",
            rw = [0.0, 0.0, 0.03],
            expected_r = 0.06,
            winding = 3,
        ),
    ]

    resulting_resistance_branches = Dict{Int,Any}()

    for case in cases
        @testset "$(case.name)" begin
            eng = deepcopy(data_eng)

            eng["transformer"]["tx1"]["rw"] = copy(case.rw)

            math = PMD.transform_data_model(
                eng;
                make_pu=false,
                kron_reduce=false,
                phase_project=false,
                correct_network_data=false,
            )

            branches = tx1_branches(math)
            transformers = tx1_transformers(math)

            println("\n========================================")
            println(case.name)
            println("rw = ", case.rw)

            println("\nVirtual transformers:")
            for (id, tr) in transformers
                println(
                    "  id=", id,
                    " source=", get(tr, "source_id", ""),
                    " f_bus=", get(tr, "f_bus", missing),
                    " t_bus=", get(tr, "t_bus", missing),
                )
            end

            println("\nTX1-derived branches:")
            for (id, br) in branches
                println(
                    "  id=", id,
                    " source=", get(br, "source_id", ""),
                    " f_bus=", br["f_bus"],
                    " t_bus=", br["t_bus"],
                    " R=", rmag(br),
                    " X=", xmag(br),
                )
            end

            # There should be exactly one TX1-derived branch carrying
            # nonzero winding resistance in each one-hot test.
            resistance_branches = [
                (id, br)
                for (id, br) in branches
                if rmag(br) > 1e-10
            ]

            @test length(resistance_branches) == 1

            if length(resistance_branches) == 1
                id, br = only(resistance_branches)

                @test rmag(br) ≈ case.expected_r atol=1e-10

                # Winding resistance branches should contain no
                # leakage reactance.
                @test xmag(br) ≈ 0.0 atol=1e-10

                resulting_resistance_branches[case.winding] = (
                    id = id,
                    source_id = get(br, "source_id", ""),
                    f_bus = br["f_bus"],
                    t_bus = br["t_bus"],
                    r = rmag(br),
                )
            end
        end
    end


    # ------------------------------------------------------------------
    # Cross-case checks
    #
    # Changing which winding has resistance should move the nonzero
    # resistance to a different mathematical branch/path.
    # ------------------------------------------------------------------
    @testset "leakage branches are correct" begin
        eng = deepcopy(data_eng)

        eng["transformer"]["tx1"]["rw"] = [0.01, 0.02, 0.03]
        eng["transformer"]["tx1"]["xsc"] = [0.04, 0.06, 0.08]

        math = PMD.transform_data_model(
            eng;
            make_pu=false,
            kron_reduce=false,
            phase_project=false,
            correct_network_data=false,
        )

        branches = tx1_branches(math)

        x_by_source = Dict(
            get(br, "source_id", "") => xmag(br)
            for (_, br) in branches
            if xmag(br) > 1e-10
        )

        @test x_by_source["_virtual_branch.transformer.tx1_4"] ≈ 0.092 atol=1e-10
        @test x_by_source["_virtual_branch.transformer.tx1_5"] ≈ (0.46/3) atol=1e-10
        @test x_by_source["_virtual_branch.transformer.tx1_6"] ≈ 0.46 atol=1e-10
    end

    
    @testset "winding identity is preserved" begin
        @test length(resulting_resistance_branches) == 3

        if length(resulting_resistance_branches) == 3
            b1 = resulting_resistance_branches[1]
            b2 = resulting_resistance_branches[2]
            b3 = resulting_resistance_branches[3]

            @test b1.source_id == "_virtual_branch.transformer.tx1_1"
            @test b2.source_id == "_virtual_branch.transformer.tx1_2"
            @test b3.source_id == "_virtual_branch.transformer.tx1_3"

            @test b1.r ≈ 0.02 atol=1e-10
            @test b2.r ≈ 0.04 atol=1e-10
            @test b3.r ≈ 0.06 atol=1e-10
        end
    end


    # ------------------------------------------------------------------
    # All resistances enabled simultaneously.
    #
    # This checks that _build_loss_model! does not lose or permute one
    # of the winding resistance branches when all three are present.
    # ------------------------------------------------------------------

    @testset "all winding resistances together" begin
        eng = deepcopy(data_eng)

        eng["transformer"]["tx1"]["rw"] = [0.01, 0.02, 0.03]

        math = PMD.transform_data_model(
            eng;
            make_pu=false,
            kron_reduce=false,
            phase_project=false,
            correct_network_data=false,
        )

        branches = tx1_branches(math)

        resistance_branches = Dict(
            get(br, "source_id", "") => rmag(br)
            for (_, br) in branches
            if rmag(br) > 1e-10
        )

        @test length(resistance_branches) == 3

        @test resistance_branches[
            "_virtual_branch.transformer.tx1_1"
        ] ≈ 0.02 atol=1e-10

        @test resistance_branches[
            "_virtual_branch.transformer.tx1_2"
        ] ≈ 0.04 atol=1e-10

        @test resistance_branches[
            "_virtual_branch.transformer.tx1_3"
        ] ≈ 0.06 atol=1e-10
        
    end
end