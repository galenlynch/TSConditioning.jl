"Apply filter b to signal x using mmaped array"
function filtfilt_mmap(
    b::AbstractVector,
    x::AbstractVector,
    basedir::AbstractString = tempdir(),
    autoclean::Bool = true
)
    nb = length(b)

    nx = length(x)
    n_exp = nb - 1
    n_offset = 2 * n_exp

    # Convolve b with its reverse, see DSP.jl filtfilt
    newb = filt(b, reverse(b))
    resize!(newb, n_offset + 1)
    for i = 1:nb-1
        newb[nb+i] = newb[nb-i]
    end

    T = Base.promote_eltype(b, x)

    # Extrapolate signal
    extrapolated = extrapolate_signal(x, n_exp)

    # Make output
    (filtered, filt_path) = typemmap(Vector{T}, (length(x) + n_offset,), basedir)

    if autoclean
        fp_copy = deepcopy(filt_path)
        atexit(() -> rm(fp_copy))
    end

    # Filter
    filt!(filtered, newb, extrapolated)

    out = @view filtered[(1 + n_offset):end]
    return (out, filt_path, n_offset)
end

"hpf signal x using mmaped array"
function filtfilt_mmap(
    x::AbstractVector,
    fs::Number,
    args...;
    fc::Real = 800.0,
    win::AbstractVector{<:AbstractFloat} = blackman(91)
)
    return filtfilt_mmap(make_hpf_taps(fc, fs, win), x, args...)
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

function hpf(
    s::AbstractArray{<:Number},
    fs::Number;
    fc::Number = 800,
    win::AbstractArray{<:AbstractFloat} = blackman(91)
)
    df = make_hpf_taps(fc, fs, win)
    return filtfilt(df, s)
end

function make_hpf_taps(
    fc::AbstractFloat, fs::Number, win::AbstractArray{<:AbstractFloat} = blackman(91)
)
    resp = Highpass(fc; fs = fs)::DSP.Filters.Highpass{Float64}
    designmethod = FIRWindow(win)
    return digitalfilter(resp, designmethod)::Vector{Float64}
end
function make_hpf_taps(fc::Number, fs, args...)
    make_hpf_taps(convert(Float64, fc), fs, args...)
end

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

function extrapolate_signal(sig::AbstractVector{E}, pad_length::Integer) where E
    if length(sig) < pad_length + 1
        throw(ArgumentError("sig length must be at least $(pad_length + 1)"))
    end
    left_pad = Vector{E}(pad_length)
    right_pad = Vector{E}(pad_length)
    xb = 2 * sig[1]
    xe = 2 * sig[end]
    @inbounds for i in 1:pad_length
        left_pad[i] = xb - sig[2 + pad_length - i]
        right_pad[i] = xe - sig[end - i]
    end
    JoinedVectors(left_pad, sig, right_pad)
end
