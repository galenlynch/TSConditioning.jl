__precompile__()
module TSConditioning

using DSP

export
    hpf,
    make_hpf_taps,
    mua,
    smooth,
    gaussian_kernel,
    rescale,
    filtfilt_mmap,
    rev_view

rev_view(a::AbstractVector) = @view a[end:-1:1]

function typemmap(::Type{A}, dims::NTuple{N, Int}) where {A<:AbstractArray, N}
    (path, io) = mktemp()
    arr = try
        Mmap.mmap(io, A, dims; grow = true)
    finally
        close(io) # Can still write to array
        rm(path) # File will remain accessible to this process
    end
    return arr
end
typemmap(a::A) where A<:AbstractArray = typemmap(A, size(a))


function filtfilt_mmap(b::AbstractVector, x::AbstractVector)
    nb = length(b)
    # Only need as much padding as the order of the filter
    # Convolve b with its reverse
    newb = filt(b, reverse(b))
    resize!(newb, 2 * nb - 1)
    for i = 1:nb-1
        newb[nb+i] = newb[nb-i]
    end

    # Extrapolate signal
    T = Base.promote_eltype(b, x)
    extrapolated = typemmap(Vector{T}, (length(x) + 2 * nb - 2,))
    DSP.Filters.extrapolate_signal!(extrapolated, 1, x, 1, length(x), nb - 1)

    # Filter
    filtered = typemmap(Vector{T}, (length(x) + 2 * nb - 2,))
    filt!(filtered, newb, extrapolated)
    out = @view filtered[(2 * nb -1):end]
    return out
end
function filtfilt_mmap(
    x::AbstractVector,
    fs::Number;
    fc::AbstractFloat = 800.0,
    win = hamming(91)
)
    return filtfilt_mmap(make_hpf_taps(fc, fs; win = win), x)
end

function hpf(s::AbstractArray, fs::Number; fc::AbstractFloat = 800, win = hamming(91))
    df = make_hpf_taps(fc, fs, win)
    return filtfilt(df, s)
end

function make_hpf_taps(fc::AbstractFloat, fs; win = hamming(91))
    resp = Highpass(fc; fs = fs)
    designmethod = FIRWindow(win)
    return digitalfilter(resp, designmethod)
end
make_hpf_taps(fc, args...; kwargs...) = make_hpf_taps(convert(Float64, fc), args...;kwargs...)

function same_conv_indices(a::Integer, b::Integer)
    offset = floor(Int, (b - 1) / 2)
    return (offset,  ndx_offset(offset, a))
end
function same_conv_indices(a::AbstractVector, b::AbstractVector)
    return same_conv_indices(length(a), length(b))
end

function gaussian_kernel(l::Integer, sig)
    basis = 0:(l - 1)
    mu = (l - 1) / 2
    k = exp.(-(basis - mu) .^ 2 / (2 * sig ^ 2))
    return k ./ sum(k)
end

function smooth(
    a::A, k::B; portion::Symbol = :same
) where {T<:Number, S<:Number, A<:AbstractVector{T}, B<:AbstractVector{S}}
    c = conv(a,k)
    if portion == :same
        (ib, ie) = same_conv_indices(a, k)
        out = c[ib:ie]
    elseif portion == :all
        out = c
    else
        error("uncrecognize portion")
    end
    return out
end
function smooth(
    a::A, k::AbstractVector; kwargs...
) where {T<:AbstractVector, A<:AbstractVector{T}}
    return map((x) -> smooth(x, k), a)
end
function smooth(a, sig::Number, fs = one(sig), l = 6*sig; kwargs...)
    sig = sig * fs
    l = ceil(Int, l * fs)
    l = isodd(l) ? l : l + 1
    k = gaussian_kernel(l, sig)
    smooth(a, k; kwargs...)
end

conv_length(inds::NTuple{2, T}) where T<: Number = sum(length.(inds)) - 1

function mua(
    a::A, sig::Number, fs::Number = 1,
    frect::Function = (x) -> x^2, args...
    ; kwargs...
) where {T<:Number, A<:AbstractVector{T}}
    rect = frect.(a)
    return smooth(rect, sig, fs, args...; kwargs...)
end
function mua(
    as::A, sig::Number, fs::Number = 1,
    frect::Function = (x) -> x.^2, args...; kwargs...
) where {T<:AbstractVector, A<:AbstractVector{T}}
    rect = similar(as)
    for (i, a) in enumerate(as)
        rect[i] = frect.(a)
    end
    return smooth(rect, sig, fs, args...; kwargs...)
end

function rescale(x, c = 1, o = 0)
    (s, b) = extrema(x)
    return c * (x - s) / (b - s) + o
end


end # module
