using Test
using MultiScales

@testset "arithmetic.jl" begin
    @testset "_halfadder" begin
        adder = MultiScales._adder(half=true)
        @test size(adder) == (1, 2, 2, 2, 2)
        for (a, b, out, cout) in [
                (0, 0, 0, 0),
                (0, 1, 1, 0),
                (1, 0, 1, 0),
                (1, 1, 0, 1),
            ]
            @test adder[1, cout+1, a+1, b+1, out+1] == 1.0
        end
    end

    @testset "_adder" begin
        adder = MultiScales._adder()
        @test size(adder) == (2, 2, 2, 2, 2)
        for (cin, a, b, out, cout) in [
                (0, 0, 0, 0, 0),
                (0, 0, 1, 1, 0),
                (0, 1, 0, 1, 0),
                (0, 1, 1, 0, 1),
                (1, 0, 0, 1, 0),
                (1, 0, 1, 0, 1),
                (1, 1, 0, 0, 1),
                (1, 1, 1, 1, 1),
            ]
            @test adder[cin + 1, cout+1, a+1, b+1, out+1] == 1.0
        end
    end

    @testset "_lastadder" begin
        for antiperiodic in [false, true]
            adder = MultiScales._adder(last=true, antiperiodic=antiperiodic)
            @test size(adder) == (2, 1, 2, 2, 2)
            for (cin, a, b, out, cout) in [
                    (0, 0, 0, 0, 0),
                    (0, 0, 1, 1, 0),
                    (0, 1, 0, 1, 0),
                    (0, 1, 1, 0, 1),
                    (1, 0, 0, 1, 0),
                    (1, 0, 1, 0, 1),
                    (1, 1, 0, 0, 1),
                    (1, 1, 1, 1, 1),
                ]
                coeff::ComplexF64 = ((antiperiodic && cout != 0) ? -1 : 1)
                @test adder[cin + 1, 1, a+1, b+1, out+1] == coeff
            end
        end
    end

    @testset "_sector" begin
        for which in [:left, :right], b in 0:1, a in 0:1
            selector = MultiScales._selector(which)
            @test size(selector) == (1, 1, 2, 2, 2)
            out = which == :left ? a : b
            @test selector[1, 1, a+1, b+1, out+1] == 1.0
        end
    end
end
