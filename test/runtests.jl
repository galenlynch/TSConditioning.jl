using Compat, TSConditioning

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test, Statistics
end

const A = rand(200)

# write your own tests here
@testset "TSConditioning" begin
    @testset "util" begin
        @test isapprox(center(A), center(A .+ 5))
        @test rescale([0, 1]) == [0, 1]
        @test rescale([1, 2]) == [0, 1]
        @test rescale([0, 1], 2) == [0, 2]
        @test rescale([1, 2], 2, 1) == [1, 3]
    end

    @testset "whitening" begin
        x = rand(3, 200)
        c = Compat.Statistics.cov(x; dims = 2)
        w = zca(c)
        y = w * (x .- Compat.Statistics.mean(x; dims = 2))
        wc = Compat.Statistics.cov(y; dims = 2)

        x1 = x[1, :]
        x2 = x[2, :]
        x3 = x[3, :]

        outs, paths = whiten_mmap([x1, x2, x3], w)
    end

    @testset "filt" begin
        make_hpf_taps(800, 2000)
        hpf(A, 2000)
    end

end
