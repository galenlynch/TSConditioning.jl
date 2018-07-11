__precompile__()
module TSConditioning

using DSP, GLUtilities, JoinedArrays

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
    xcorr_centered,
    norm_sig_xcorr!,
    norm_sig_xcorr,
    whiten_mmap,
    xcorr_normed,
    xcorr_ndx_lag,
    xcorr_lags,
    xcorr_best_lag,
    zcd

include("util.jl")
include("filt.jl")
include("xcorr.jl")
include("whitening.jl")
end # module
