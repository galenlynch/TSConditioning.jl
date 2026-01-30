__precompile__()
module TSConditioning

using DSP: DSP, Bandpass, FIRFilter, FIRWindow, Highpass, blackman, conv,
    digitalfilter, filt, filt!, filtfilt, xcorr
using GLUtilities: GLUtilities, allsame, bin_bounds, div_type, n_ndx, ndx_offset,
    rev_view, typemmap
using JoinedArrays: JoinedArrays, JoinedVectors
using LinearAlgebra: LinearAlgebra, Diagonal, Symmetric, mul!, svd
using MultivariateStats: MultivariateStats, llsq
using Statistics: Statistics, cov, mean, std

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
