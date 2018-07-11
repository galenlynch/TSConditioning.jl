center(x::AbstractArray) = x .- mean(x)

"Remove minimum value, scale, and move signal"
function rescale(x, c = 1, o = 0)
    (s, b) = extrema(x)
    return c * (x - s) / (b - s) + o
end
