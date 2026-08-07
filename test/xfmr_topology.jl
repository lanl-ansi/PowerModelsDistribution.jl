using Test
using PowerModelsDistribution

const PMD = PowerModelsDistribution

@testset "transformer topology and winding association" begin

    file = "../test/data/opendss/ut_trans_3w_dyy_basetest.dss"
    data_eng = open(file) do io
       parse_opendss(io)
       end
       
    data_eng["transformer"]["tx1"]["xsc"] = [0.04, 0.06, 0.08]


    function build_math(eng)
        PMD.transform_data_model(
            eng;
            make_pu=false,
            kron_reduce=false,
            phase_project=false,
            correct_network_data=false,
        )
    end

    function tx1_transformers(math)
        Dict(
            get(tr, "source_id", "") => tr
            for (_, tr) in math["transformer"]
            if occursin(
                "transformer.tx1",
                lowercase(string(get(tr, "source_id", "")))
            )
        )
    end

    function tx1_branches(math)
        Dict(
            get(br, "source_id", "") => br
            for (_, br) in math["branch"]
            if occursin(
                "transformer.tx1",
                lowercase(string(get(br, "source_id", "")))
            )
        )
    end

    rmag(br) = maximum(abs, br["br_r"])
    xmag(br) = maximum(abs, br["br_x"])

    # Does a branch touch a given mathematical bus?
    touches(br, bus) = br["f_bus"] == bus || br["t_bus"] == bus

    # Return the bus on the opposite side of a branch.
    function other_bus(br, bus)
        if br["f_bus"] == bus
            return br["t_bus"]
        elseif br["t_bus"] == bus
            return br["f_bus"]
        else
            error("branch does not touch bus $bus")
        end
    end


    @testset "each winding maps to its own resistance path" begin
        eng = deepcopy(data_eng)

        tx = eng["transformer"]["tx1"]

        tx["rw"] = [0.01, 0.02, 0.03]
        tx["xsc"] = [0.04, 0.06, 0.08]

        math = build_math(eng)

        trs = tx1_transformers(math)
        brs = tx1_branches(math)

        # Expected derived object names.
        tr1_id = "_virtual_transformer.transformer.tx1.1"
        tr2_id = "_virtual_transformer.transformer.tx1.2"
        tr3_id = "_virtual_transformer.transformer.tx1.3"

        r1_id = "_virtual_branch.transformer.tx1_1"
        r2_id = "_virtual_branch.transformer.tx1_2"
        r3_id = "_virtual_branch.transformer.tx1_3"

        @test haskey(trs, tr1_id)
        @test haskey(trs, tr2_id)
        @test haskey(trs, tr3_id)

        @test haskey(brs, r1_id)
        @test haskey(brs, r2_id)
        @test haskey(brs, r3_id)

        tr1 = trs[tr1_id]
        tr2 = trs[tr2_id]
        tr3 = trs[tr3_id]

        r1 = brs[r1_id]
        r2 = brs[r2_id]
        r3 = brs[r3_id]

        # Numerical checks.
        @test rmag(r1) ≈ 0.02 atol=1e-10
        @test rmag(r2) ≈ 0.04 atol=1e-10
        @test rmag(r3) ≈ 0.06 atol=1e-10

        @test xmag(r1) ≈ 0.0 atol=1e-10
        @test xmag(r2) ≈ 0.0 atol=1e-10
        @test xmag(r3) ≈ 0.0 atol=1e-10

        # The winding resistance branch must touch the internal
        # t_bus of the matching virtual transformer.
        @test touches(r1, tr1["t_bus"])
        @test touches(r2, tr2["t_bus"])
        @test touches(r3, tr3["t_bus"])

        # It should NOT touch the internal winding-side bus of another
        # winding.
        @test !touches(r1, tr2["t_bus"])
        @test !touches(r1, tr3["t_bus"])

        @test !touches(r2, tr1["t_bus"])
        @test !touches(r2, tr3["t_bus"])

        @test !touches(r3, tr1["t_bus"])
        @test !touches(r3, tr2["t_bus"])

        println("\nWinding-to-resistance topology:")
        for (n, tr, br) in [
            (1, tr1, r1),
            (2, tr2, r2),
            (3, tr3, r3),
        ]
            println(
                "winding ", n,
                " original f_bus=", tr["f_bus"],
                " internal bus=", tr["t_bus"],
                " resistance branch=(",
                br["f_bus"], ",", br["t_bus"], ")",
                " R=", rmag(br),
            )
        end
    end


    @testset "leakage network joins the three winding resistance paths" begin
        eng = deepcopy(data_eng)

        tx = eng["transformer"]["tx1"]

        tx["rw"] = [0.01, 0.02, 0.03]
        tx["xsc"] = [0.04, 0.06, 0.08]

        math = build_math(eng)

        trs = tx1_transformers(math)
        brs = tx1_branches(math)

        r1 = brs["_virtual_branch.transformer.tx1_1"]
        r2 = brs["_virtual_branch.transformer.tx1_2"]
        r3 = brs["_virtual_branch.transformer.tx1_3"]

        # Find the leakage-side node of each resistance branch.
        #
        # One endpoint is the matching virtual-transformer t_bus.
        # The other endpoint belongs to the internal leakage network.
        leakage_node_1 = other_bus(
            r1,
            trs["_virtual_transformer.transformer.tx1.1"]["t_bus"],
        )

        leakage_node_2 = other_bus(
            r2,
            trs["_virtual_transformer.transformer.tx1.2"]["t_bus"],
        )

        leakage_node_3 = other_bus(
            r3,
            trs["_virtual_transformer.transformer.tx1.3"]["t_bus"],
        )

        @test leakage_node_1 != leakage_node_2
        @test leakage_node_1 != leakage_node_3
        @test leakage_node_2 != leakage_node_3

        # Expected leakage branch identities for a 3-winding transformer.
        z12 = brs["_virtual_branch.transformer.tx1_4"]
        z13 = brs["_virtual_branch.transformer.tx1_5"]
        z23 = brs["_virtual_branch.transformer.tx1_6"]

        # These should connect the three leakage-side winding nodes
        # pairwise.
        @test Set([z12["f_bus"], z12["t_bus"]]) ==
              Set([leakage_node_1, leakage_node_2])

        @test Set([z13["f_bus"], z13["t_bus"]]) ==
              Set([leakage_node_1, leakage_node_3])

        @test Set([z23["f_bus"], z23["t_bus"]]) ==
              Set([leakage_node_2, leakage_node_3])

        # Numerical leakage values from xsc=[0.04,0.06,0.08].
        @test xmag(z12) ≈ 0.092 atol=1e-10
        @test xmag(z13) ≈ (0.46 / 3) atol=1e-10
        @test xmag(z23) ≈ 0.46 atol=1e-10

        println("\nLeakage-side winding nodes:")
        println("w1 -> ", leakage_node_1)
        println("w2 -> ", leakage_node_2)
        println("w3 -> ", leakage_node_3)

        println("\nLeakage branches:")
        println("1-2: ", z12["f_bus"], " <-> ", z12["t_bus"])
        println("1-3: ", z13["f_bus"], " <-> ", z13["t_bus"])
        println("2-3: ", z23["f_bus"], " <-> ", z23["t_bus"])
    end


    @testset "winding permutation preserves topology" begin
        # --------------------------------------------------------------
        # Original:
        #
        # winding 1: bus 1, 11 kV, rw=0.01
        # winding 2: bus 2,  4 kV, rw=0.02
        # winding 3: bus 3, .4 kV, rw=0.03
        #
        # x12=0.04, x13=0.06, x23=0.08
        #
        # Now swap winding labels 2 <-> 3 consistently.
        # --------------------------------------------------------------

        eng_a = deepcopy(data_eng)

        txa = eng_a["transformer"]["tx1"]
        txa["rw"] = [0.01, 0.02, 0.03]
        txa["xsc"] = [0.04, 0.06, 0.08]

        math_a = build_math(eng_a)

        eng_b = deepcopy(data_eng)
        txb = eng_b["transformer"]["tx1"]

        # Swap winding 2 and 3 in every winding-indexed field.
        perm = [1, 3, 2]

        for key in [
            "bus",
            "connections",
            "configuration",
            "polarity",
            "vm_nom",
            "sm_nom",
            "tm_set",
            "tm_lb",
            "tm_ub",
            "tm_fix",
            "tm_step",
        ]
            txb[key] = txb[key][perm]
        end

        txb["rw"] = [0.01, 0.03, 0.02]

        # Original:
        #   x12 = 0.04
        #   x13 = 0.06
        #   x23 = 0.08
        #
        # After swapping winding labels 2 <-> 3:
        #   new x12 = old x13 = 0.06
        #   new x13 = old x12 = 0.04
        #   new x23 = old x23 = 0.08
        txb["xsc"] = [0.06, 0.04, 0.08]

        math_b = build_math(eng_b)

        brs_a = tx1_branches(math_a)
        brs_b = tx1_branches(math_b)

        # Original resistance branch values.
        @test rmag(brs_a["_virtual_branch.transformer.tx1_1"]) ≈ 0.02 atol=1e-10
        @test rmag(brs_a["_virtual_branch.transformer.tx1_2"]) ≈ 0.04 atol=1e-10
        @test rmag(brs_a["_virtual_branch.transformer.tx1_3"]) ≈ 0.06 atol=1e-10

        # After permutation, winding labels 2 and 3 have swapped,
        # so the derived source IDs should carry the correspondingly
        # swapped resistance values.
        @test rmag(brs_b["_virtual_branch.transformer.tx1_1"]) ≈ 0.02 atol=1e-10
        @test rmag(brs_b["_virtual_branch.transformer.tx1_2"]) ≈ 0.06 atol=1e-10
        @test rmag(brs_b["_virtual_branch.transformer.tx1_3"]) ≈ 0.04 atol=1e-10

        # Leakage pair identities should swap in the same way.
        #
        # A:
        #   branch 4 = pair 1-2
        #   branch 5 = pair 1-3
        #   branch 6 = pair 2-3
        #
        # B:
        #   pair 1-2 corresponds physically to old 1-3
        #   pair 1-3 corresponds physically to old 1-2
        #   pair 2-3 is unchanged physically.

        @test xmag(brs_b["_virtual_branch.transformer.tx1_4"]) ≈
              xmag(brs_a["_virtual_branch.transformer.tx1_5"]) atol=1e-10

        @test xmag(brs_b["_virtual_branch.transformer.tx1_5"]) ≈
              xmag(brs_a["_virtual_branch.transformer.tx1_4"]) atol=1e-10

        @test xmag(brs_b["_virtual_branch.transformer.tx1_6"]) ≈
              xmag(brs_a["_virtual_branch.transformer.tx1_6"]) atol=1e-10
    end

end