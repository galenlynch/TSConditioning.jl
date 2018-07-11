function zca(c::AbstractArray{T, 2}) where T
    R = div_type(T)
    F = svdfact(c)
    u = F[:U]::Matrix{R}
    s = F[:S]::Vector{R}
    u * diagm(s .^ -(1 / 2)) * u'
end

function whiten_mmap(
    xs::AbstractVector{<:AbstractVector{T}}, w::AbstractArray{S, 2}
) where {T, S}
    nx = length(xs)
    nx > 0 || throw(ArgumentError("xs is empty"))
    nel = length(xs[1])
    allsame(length, xs) || throw(ArgumentError("xs not the same length"))
    w_sz = size(w)
    w_sz == (nx, nx) || throw(
        ArgumentError("Size of covariance matrix does not match")
    )
    R = promote_type(T, S)
    outs = Vector{Vector{R}}(nx)
    paths = Vector{String}(nx)
    for i = 1:nx
        outs[i], paths[i] = typemmap(Vector{R}, (nel,))
    end
    _whiten_mmap!(outs, xs, w, nx, nel)
    outs, paths
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
