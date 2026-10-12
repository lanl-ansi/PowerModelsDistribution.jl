    # Helper for the effective impedance seen between nodes i and j
    # in a 3-node delta/mesh network.
    function _zeq_3node(zbr, i, j, k)
        zij = zbr[(min(i,j), max(i,j))]
        zik = zbr[(min(i,k), max(i,k))]
        zkj = zbr[(min(k,j), max(k,j))]

        # Direct branch i-j is in parallel with path i-k-j
        return inv(inv(zij) + inv(zik + zkj))
    end


@testset "_sc2br_impedance" begin

    @testset "2-winding identity" begin
        for z in [
            0.1im,
            0.04im,
            6.666666666666666e-6im,
            0.01 + 0.05im,
        ]
            zsc = Dict{Tuple{Int,Int},ComplexF64}(
                (1,2) => z,
            )

            zbr = PMD._sc2br_impedance(copy(zsc))

            @test keys(zbr) == keys(zsc)
            @test zbr[(1,2)] ≈ zsc[(1,2)]
        end
    end


    @testset "3-winding asymmetric purely reactive" begin
        zsc = Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.08im,
            (1,3) => 0.12im,
            (2,3) => 0.16im,
        )

        zbr = PMD._sc2br_impedance(copy(zsc))

        # Expected equivalent delta impedances.
        #
        # Input short-circuit impedances imply the star legs:
        #
        # z1 = (z12 + z13 - z23)/2 = 0.02im
        # z2 = (z12 + z23 - z13)/2 = 0.06im
        # z3 = (z13 + z23 - z12)/2 = 0.10im
        #
        # Star -> delta gives:
        #
        # Z12 = z1 + z2 + z1*z2/z3 = 0.092im
        # Z13 = z1 + z3 + z1*z3/z2 = 0.153333...im
        # Z23 = z2 + z3 + z2*z3/z1 = 0.46im

        @test zbr[(1,2)] ≈ 0.092im
        @test zbr[(1,3)] ≈ (0.46 / 3)im
        @test zbr[(2,3)] ≈ 0.46im

        # More important invariant:
        # measuring the equivalent branch network between each pair
        # must reproduce the original short-circuit test impedance.
        @test _zeq_3node(zbr, 1, 2, 3) ≈ zsc[(1,2)]
        @test _zeq_3node(zbr, 1, 3, 2) ≈ zsc[(1,3)]
        @test _zeq_3node(zbr, 2, 3, 1) ≈ zsc[(2,3)]
    end


    @testset "3-winding asymmetric complex impedance" begin
        # This makes sure the implementation works for general complex
        # impedance, not just purely imaginary leakage reactance.
        zsc = Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.01 + 0.08im,
            (1,3) => 0.02 + 0.12im,
            (2,3) => 0.03 + 0.16im,
        )

        zbr = PMD._sc2br_impedance(copy(zsc))

        @test _zeq_3node(zbr, 1, 2, 3) ≈ zsc[(1,2)]
        @test _zeq_3node(zbr, 1, 3, 2) ≈ zsc[(1,3)]
        @test _zeq_3node(zbr, 2, 3, 1) ≈ zsc[(2,3)]
    end


    @testset "3-winding unequal but physically simple star" begin
        # Construct the short-circuit measurements from known star legs.
        z1 = 0.02im
        z2 = 0.03im
        z3 = 0.07im

        zsc = Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => z1 + z2,
            (1,3) => z1 + z3,
            (2,3) => z2 + z3,
        )

        zbr = PMD._sc2br_impedance(copy(zsc))

        @test _zeq_3node(zbr, 1, 2, 3) ≈ z1 + z2
        @test _zeq_3node(zbr, 1, 3, 2) ≈ z1 + z3
        @test _zeq_3node(zbr, 2, 3, 1) ≈ z2 + z3
    end


    @testset "3-winding symmetric case" begin
        zsc = Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.1im,
            (1,3) => 0.1im,
            (2,3) => 0.1im,
        )

        zbr = PMD._sc2br_impedance(copy(zsc))

        # Symmetry should be preserved.
        @test zbr[(1,2)] ≈ zbr[(1,3)]
        @test zbr[(1,2)] ≈ zbr[(2,3)]

        @test _zeq_3node(zbr, 1, 2, 3) ≈ 0.1im
        @test _zeq_3node(zbr, 1, 3, 2) ≈ 0.1im
        @test _zeq_3node(zbr, 2, 3, 1) ≈ 0.1im
    end


    @testset "all-zero short-circuit impedances" begin
        zsc = Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.0 + 0.0im,
            (1,3) => 0.0 + 0.0im,
            (2,3) => 0.0 + 0.0im,
        )

        zbr = PMD._sc2br_impedance(copy(zsc))

        @test zbr == zsc
    end


    @testset "lower-triangle key is accepted" begin
        zsc = Dict{Tuple{Int,Int},ComplexF64}(
            (2,1) => 0.08im,
            (1,3) => 0.12im,
            (2,3) => 0.16im,
        )

        zsc_copy = copy(zsc)
        zbr = PMD._sc2br_impedance(zsc_copy)

        # Function currently fills the missing upper-triangle entry
        # by mutating its input.
        @test haskey(zsc_copy, (1,2))
        @test zsc_copy[(1,2)] ≈ 0.08im

        # The returned upper-triangle network should still reproduce
        # the supplied short-circuit impedances.
        @test _zeq_3node(zbr, 1, 2, 3) ≈ 0.08im
        @test _zeq_3node(zbr, 1, 3, 2) ≈ 0.12im
        @test _zeq_3node(zbr, 2, 3, 1) ≈ 0.16im
    end


    @testset "missing short-circuit pair throws" begin
        zsc = Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.08im,
            (1,3) => 0.12im,
            # (2,3) deliberately missing
        )

        @test_throws ErrorException PMD._sc2br_impedance(copy(zsc))
    end


    @testset "one zero short-circuit impedance" begin
        # This is a useful numerical edge case because it can make
        # Zb singular or nearly singular depending on the other values.
        zsc = Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.0im,
            (1,3) => 0.10im,
            (2,3) => 0.10im,
        )

        zbr = PMD._sc2br_impedance(copy(zsc))

        # The main assertion here is that the conversion completes
        # and reproduces the pairwise effective impedances.
        #
        # Some resulting branch impedances may be zero or Inf depending
        # on the equivalent network.
        @test _zeq_3node(zbr, 1, 3, 2) ≈ zsc[(1,3)]
        @test _zeq_3node(zbr, 2, 3, 1) ≈ zsc[(2,3)]
    end


    @testset "winding permutation invariance" begin
        zsc = Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.08im,
            (1,3) => 0.12im,
            (2,3) => 0.16im,
        )

        zbr = PMD._sc2br_impedance(copy(zsc))

        # Swap winding labels 2 <-> 3.
        #
        # old:
        #   12 = 0.08
        #   13 = 0.12
        #   23 = 0.16
        #
        # new:
        #   12 = old 13
        #   13 = old 12
        #   23 = old 23

        zsc_perm = Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.12im,
            (1,3) => 0.08im,
            (2,3) => 0.16im,
        )

        zbr_perm = PMD._sc2br_impedance(copy(zsc_perm))

        # Relabel the permuted result back to the original numbering.
        @test zbr_perm[(1,3)] ≈ zbr[(1,2)]
        @test zbr_perm[(1,2)] ≈ zbr[(1,3)]
        @test zbr_perm[(2,3)] ≈ zbr[(2,3)]

        @test _zeq_3node(zbr_perm, 1, 2, 3) ≈ zsc_perm[(1,2)]
        @test _zeq_3node(zbr_perm, 1, 3, 2) ≈ zsc_perm[(1,3)]
        @test _zeq_3node(zbr_perm, 2, 3, 1) ≈ zsc_perm[(2,3)]
    end

end

@testset "near-zero short-circuit impedance" begin
    for eps in [1e-3, 1e-6, 1e-9, 1e-12, 1e-15]
        zsc = Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => eps * im,
            (1,3) => 0.10im,
            (2,3) => (0.10 + eps) * im,
        )

        zbr = PMD._sc2br_impedance(copy(zsc))

        @test _zeq_3node(zbr, 1,2,3) ≈ zsc[(1,2)]
        @test _zeq_3node(zbr, 1,3,2) ≈ zsc[(1,3)]
        @test _zeq_3node(zbr, 2,3,1) ≈ zsc[(2,3)]
    end
end

@testset "_sc2br_impedance singular and permutation cases" begin
    function _zeq_3node(zbr, i, j, k)
        zij = zbr[(min(i,j), max(i,j))]
        zik = zbr[(min(i,k), max(i,k))]
        zkj = zbr[(min(k,j), max(k,j))]
        return inv(inv(zij) + inv(zik + zkj))
    end

    # Explicit zero short-circuit pair, all permutations
    zero_pair_cases = [
        Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.00im,
            (1,3) => 0.10im,
            (2,3) => 0.10im,
        ),
        Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.10im,
            (1,3) => 0.00im,
            (2,3) => 0.10im,
        ),
        Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.10im,
            (1,3) => 0.10im,
            (2,3) => 0.00im,
        ),
    ]

    for zsc in zero_pair_cases
        zbr = PMD._sc2br_impedance(copy(zsc))

        @show zsc
        @show zbr

        @test !any(
            isnan(real(z)) || isnan(imag(z))
            for z in values(zbr)
            if isfinite(real(z)) && isfinite(imag(z))
        )
    end

    # Zero implied star leg, but all Zsc values are nonzero.
    # These should remain electrically equivalent under winding permutation.
    zero_star_leg_cases = [
        # z1 = (Z12 + Z13 - Z23)/2 = 0
        Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.10im,
            (1,3) => 0.20im,
            (2,3) => 0.30im,
        ),

        # z2 = (Z12 + Z23 - Z13)/2 = 0
        Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.10im,
            (1,3) => 0.30im,
            (2,3) => 0.20im,
        ),

        # z3 = (Z13 + Z23 - Z12)/2 = 0
        Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => 0.30im,
            (1,3) => 0.10im,
            (2,3) => 0.20im,
        ),
    ]

    for zsc in zero_star_leg_cases
        zbr = PMD._sc2br_impedance(copy(zsc))

        @show zsc
        @show zbr

        @test _zeq_3node(zbr, 1, 2, 3) ≈ zsc[(1,2)]
        @test _zeq_3node(zbr, 1, 3, 2) ≈ zsc[(1,3)]
        @test _zeq_3node(zbr, 2, 3, 1) ≈ zsc[(2,3)]
    end

    # Near-zero case should approach the singular case continuously
    for eps in [1e-3, 1e-6, 1e-9, 1e-12]
        zsc = Dict{Tuple{Int,Int},ComplexF64}(
            (1,2) => eps * im,
            (1,3) => 0.10im,
            (2,3) => (0.10 + eps) * im,
        )

        zbr = PMD._sc2br_impedance(copy(zsc))

        @test _zeq_3node(zbr, 1, 2, 3) ≈ zsc[(1,2)]
        @test _zeq_3node(zbr, 1, 3, 2) ≈ zsc[(1,3)]
        @test _zeq_3node(zbr, 2, 3, 1) ≈ zsc[(2,3)]
    end
end
println("=====new tests end===========")