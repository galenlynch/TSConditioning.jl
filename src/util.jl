"Remove mean from signal"
center!(x::AbstractArray) = x .-= mean(x)

"See `center!`"
center(x::AbstractArray) = center!(copy(x))

"Fit linear trend to signal"
function lin_trend(ys::AbstractVector{T}) where T<:AbstractFloat
    ny = length(ys)
    basis = reshape(convert(Vector{T}, 1:ny), ny, 1)
    llsq(basis, ys)
end

"Remove linear trend from signal"
function detrend!(ys, a = lin_trend(ys))
    length(a) > 1 || throw(ArgumentError("a must have at least two elements"))
    @inbounds @simd for i = 1:length(ys)
        ys[i] -= i * a[1] + a[2]
    end
    ys
end

"See `detrend!`"
detrend(ys) = detrend!(copy(ys))

"Remove minimum value, scale, and move signal"
function rescale!(x, c = 1, o = 0)
    (s, b) = extrema(x)
    x .= c .* (x .- s) ./ (b - s) .+ o
end

"See `rescale!`"
rescale(x, args...) = rescale!(copy(x), args...)
