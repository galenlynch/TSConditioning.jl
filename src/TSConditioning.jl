__precompile__()
module TSConditioning

using Compat, DSP, GLUtilities, JoinedArrays

@static if VERSION >= v"0.7.0-DEV.2575"
    using Statistics, LinearAlgebra
end

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
    center,
    xcorr_centered,
    norm_sig_xcorr!,
    norm_sig_xcorr,
    whiten_mmap,
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
