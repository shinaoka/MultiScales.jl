using Test
using MultiScales

@testset "arithmetic.jl" begin
    @testset "_halfader" begin
        ader = MultiScales._ader(half=true)
        @test size(ader) == (1, 2, 2, 2, 2)
        for (a, b, out, cout) in [
                (0, 0, 0, 0),
                (0, 1, 1, 0),
                (1, 0, 1, 0),
                (1, 1, 0, 1),
            ]
            @test ader[1, cout+1, a+1, b+1, out+1] == 1.0
        end
    end

    @testset "_ader" begin
        ader = MultiScales._ader()
        @test size(ader) == (2, 2, 2, 2, 2)
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
            @test ader[cin + 1, cout+1, a+1, b+1, out+1] == 1.0
        end
    end

    @testset "_lastader" begin
        for antiperiodic in [false, true]
            ader = MultiScales._ader(last=true, antiperiodic=antiperiodic)
            @test size(ader) == (2, 1, 2, 2, 2)
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
                @test ader[cin + 1, 1, a+1, b+1, out+1] == coeff
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
