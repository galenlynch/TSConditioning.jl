using TSConditioning
using Base.Test

const A = rand(200)

# write your own tests here
@testset "TSConditioning" begin
    @testset "util" begin
        @test isapprox(center(A), center(A + 5))
        @test all(A[end:-1:1] .== rev_view(A))
        @test rescale([0, 1]) == [0, 1]
        @test rescale([1, 2]) == [0, 1]
        @test rescale([0, 1], 2) == [0, 2]
        @test rescale([1, 2], 2, 1) == [1, 3]
        (arr, path) = TSConditioning.typemmap(Vector{Int}, (2,); cleanup=false)
        try
            @test typeof(arr) == Vector{Int}
        finally
            rm(path)
        end
    end

    @teset "filt" begin
        make_hpf_taps(800, 2000)
        hpf(A, 2000)
    end
end
