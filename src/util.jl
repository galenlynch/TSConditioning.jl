center(x::AbstractArray) = x .- mean(x)

rev_view(a::AbstractVector) = @view a[end:-1:1]

function typemmap(::Type{A}, dims::NTuple{N, Int}, basedir::AbstractString = tempdir()) where {A<:AbstractArray, N}
    (path, io) = mktemp(basedir)
    arr = try
        Mmap.mmap(io, A, dims; grow = true)
    finally
        close(io)
    end
    return (arr, path)
end
typemmap(a::A) where A<:AbstractArray = typemmap(A, size(a))

function rescale(x, c = 1, o = 0)
    (s, b) = extrema(x)
    return c * (x - s) / (b - s) + o
end
