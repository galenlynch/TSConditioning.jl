__precompile__()
module TSConditioning

using DSP, GLUtilities

export
    hpf,
    make_hpf_taps,
    mua,
    smooth,
    gaussian_kernel,
    rescale,
    filtfilt_mmap,
    filtfilt_mmap_path,
    filtfilt_stream!,
    filtfilt_stream,
    filtfilt_stream_path,
    rev_view,
    center,
    cxcorr,
    norm_sig_xcorr!,
    norm_sig_xcorr,
    normxcorr

center(x::AbstractArray) = x .- mean(x)

cxcorr(u, v) = xcorr(center(u), center(v))

function norm_sig_xcorr!(a, u)
    um = mean(u)
    uvar = std(u; corrected = false) * sqrt(length(u))
    return (u .- um) ./ uvar
end

function norm_sig_xcorr(u)
    a = similar(u)
    return norm_sig_xcorr!(a, u)
end

normxcorr(u, v) = xcorr(norm_sig_xcorr(u), norm_sig_xcorr(v))

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

function filtfilt_mmap(b::AbstractVector, x::AbstractVector, basedir::AbstractString = tempdir())
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
    exp_len = length(x) + 2 * nb - 2;
    (extrapolated, extr_path) = typemmap(Vector{T}, (exp_len,), basedir)
    (filtered, filt_path) = try
        DSP.Filters.extrapolate_signal!(extrapolated, 1, x, 1, length(x), nb - 1)

        # Filter
        (filtered, filt_path) = typemmap(Vector{T}, (length(x) + 2 * nb - 2,), basedir)
        atexit(() -> rm(filt_path))
        filt!(filtered, newb, extrapolated)
        (filtered, filt_path)
    finally
        rm(extr_path)
    end
    offset = 2 * nb - 2
    out = @view filtered[(1 + offset):end]
    return (out, filt_path, offset)
end
function filtfilt_mmap(
    x::AbstractVector,
    fs::Number,
    args...;
    fc::AbstractFloat = 800.0,
    win = hamming(91)
)
    return filtfilt_mmap(make_hpf_taps(fc, fs; win = win), x, args...)
end

function filtfilt_mmap_path(args...; kwargs...)
    (arr, path, offset) = filtfilt_mmap(args...; kwargs...)
    l = length(arr)
    adj_l = l - offset
    return (path, typeof(arr), adj_l, offset * sizeof(eltype(arr)))
end

function filtfilt_stream!(
    sout::TOut,
    f::FIRFilter,
    sigin::TIn;
    blocksize::Integer = 1024
) where {TOut <: AbstractVector, TIn <: AbstractVector}
    # Condition filter
    nb = length(f.h)
    nh = length(f.h)
    @assert length(sout) > nh "Must be at least nh"
    buff = TOut(blocksize)
    _filtfilt_stream!(sout, f, sigin, nh, buff)
    _filtfilt_stream!(rev_view(sout), f, rev_view(sout), nh, buff)
end

function _filtfilt_stream!(sout, f, sigin, nh, buff)
    blocksize = length(buff)
    ns = length(sigin)
    # warm filter state to steady state
    filt!(view(buff, 1:nh), f, 2 * sigin[1] .+ view(sigin, (nh + 1):-1:2))

    # Filter the data in blocks
    nfullbin = fld(ns, blocksize)
    for binno in 1:nfullbin
        bbounds = bin_bounds(binno, blocksize)
        filt!(buff, f, view(sigin, bbounds[1]:bbounds[2]))
        sout[bbounds[1]:bbounds[2]] = buff
    end
    if ns % blocksize > 0
        bbounds = bin_bounds(nfullbin + 1, blocksize, ns)
        ntrail = n_ndx(bbounds...)
        filt!(view(buff, 1:ntrail), f, view(sigin, bbounds[1]:bbounds[2]))
        sout[bbounds[1]:bbounds[2]] = buff[1:ntrail]
    end
end

function filtfilt_stream(f::FIRFilter, x::AbstractVector, basedir::AbstractString = pwd(); kwargs...)
    T = Base.promote_eltype(f.h, x)
    (out, out_path) = typemmap(Vector{T}, (length(x),), basedir)
    filtfilt_stream!(out, f, x; kwargs...)
    return (out, out_path)
end
filtfilt_stream(f::AbstractVector, args...; kwargs...) = filtfilt_stream(FIRFilter(f), args...; kwargs...)

function filtfilt_stream_path(args...; kwargs...)
    (out, path) = filtfilt_stream(args...; kwargs...)
    return (path, typeof(out), length(out))
end

function hpf(s::AbstractArray, fs::Number; fc::AbstractFloat = 800, win = hamming(91))
    df = make_hpf_taps(fc, fs; win = win)
    return filtfilt(df, s)
end

function make_hpf_taps(fc::AbstractFloat, fs; win = hamming(91))
    resp = Highpass(fc; fs = fs)
    designmethod = FIRWindow(win)
    return digitalfilter(resp, designmethod)
end
make_hpf_taps(fc, fs; kwargs...) = make_hpf_taps(convert(Float64, fc), fs; kwargs...)

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
