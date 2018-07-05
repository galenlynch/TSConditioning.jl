center(x::AbstractArray) = x .- mean(x)

rev_view(a::AbstractVector) = @view a[end:-1:1]

"Make a mmaped array of type A"
function typemmap(
    ::Type{A}, dims::NTuple{N, Int}, basedir::AbstractString = tempdir();
    cleanup::Bool = true
) where {A<:AbstractArray, N}
    (path, io) = mktemp(basedir)
    cleanup && atexit(()->rm(path))
    arr = try
        Mmap.mmap(io, A, dims; grow = true)
    finally
        close(io)
    end
    return (arr, path::String)
end
typemmap(a::A; kwargs...) where A<:AbstractArray = typemmap(A, size(a); kwargs...)

"Remove minimum value, scale, and move signal"
function rescale(x, c = 1, o = 0)
    (s, b) = extrema(x)
    return c * (x - s) / (b - s) + o
end
