__precompile__()
module TSConditioning

using Compat, DSP, GLUtilities, JoinedArrays, Statistics, MultivariateStats

@static if VERSION >= v"0.7.0-DEV.2575"
    using Statistics, LinearAlgebra
end

export center!,
    center,
    detrend!,
    detrend,
    hpf,
    lin_trend,
    make_hpf_taps,
    make_bandpass,
    mua,
    norm_sig_xcorr!,
    norm_sig_xcorr,
    filtfilt_mmap,
    filtfilt_mmap_path,
    filtfilt_stream!,
    filtfilt_stream,
    filtfilt_stream_path,
    gaussian_kernel,
    rescale!,
    rescale,
    smooth,
    whiten_mmap,
    xcorr_centered,
    xcorr_normed,
    xcorr_ndx_lag,
    xcorr_lags,
    xcorr_best_lag,
    zca

include("util.jl")
include("filt.jl")
include("xcorr.jl")
include("whitening.jl")
end # module
