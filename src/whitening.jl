function zca(c::AbstractArray{T, 2}) where T
    R = div_type(T)
    F = svdfact(c)
    u = F[:U]::Matrix{R}
    s = F[:S]::Vector{R}
    convert(Matrix{R}, u * diagm(s .^ -(1 / 2)) * u')
end

function whiten_mmap(
    xs::AbstractVector{<:AbstractVector{T}},
    w::AbstractArray{S, 2},
    basedir::AbstractString = tempdir(),
    fids::AbstractVector{<:Integer} = Int[],
    autoclean::Bool = true,
    suffix::AbstractString = "_whitened"
) where {T, S}
    # Input validation
    nx = length(xs)
    nx > 0 || throw(ArgumentError("xs is empty"))
    nel = length(xs[1])
    allsame(length, xs) || throw(ArgumentError("xs not the same length"))
    w_sz = size(w)
    w_sz == (nx, nx) || throw(
        ArgumentError("Size of covariance matrix does not match")
    )
    if ! isempty(fids) && length(fids) != nx
        throw(ArgumentError("fids must be empty or same size as xs"))
    end

    # Pre-allocation of ouputs
    R = promote_type(T, S)
    outs = Vector{Vector{R}}(nx)
    paths = Vector{String}(nx)
    if isempty(fids)
        @. paths = joinpath(
            basedir,
            string(
                basename(tempname()), suffix, '_', R, '_', nel, ".dat"
            )
        )
    else
        @. paths = joinpath(basedir, string(fids, suffix, '_', R, '_', nel, ".dat"))
    end

    # Make array to be mutated
    for i = 1:nx
        outs[i], paths[i] = typemmap(
            Vector{R}, (nel,); fpath = paths[i], autoclean = autoclean
        )
    end

    # Whiten
    _whiten_mmap!(outs, xs, w, nx, nel)

    outs, paths
end

function whiten_mmap(
    xs::AbstractVector{<:AbstractVector{T}},
    basedir::AbstractString = tempdir(),
    args...
) where T
    C = collect(Symmetric(cov(xs)))
    W = zca(C)
    outs, paths = whiten_mmap(xs, W, basedir, args...)
    outs, paths, W, C
end

# assumes outs has at least one element
# Does no error checking
function _whiten_mmap!(
    outs::AbstractVector{<:AbstractVector{R}},
    xs::AbstractVector{<:AbstractVector{T}},
    w::AbstractArray{<:Number, 2},
    nx::Integer,
    nel::Integer
) where {T, R}
    means = Vector{div_type(T)}(nx)
    scratch = Vector{R}(nx)
    whitened = Vector{R}(nx)
    means .= mean.(xs)
    for elno = 1:nel
        @simd for xno = 1:nx
            @inbounds scratch[xno] = xs[xno][elno] - means[xno]
        end
        A_mul_B!(whitened, w, scratch)
        @simd for xno = 1:nx
            @inbounds outs[xno][elno] = whitened[xno]
        end
    end
    nothing
end
