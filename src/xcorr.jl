xcorr_centered(u, v) = xcorr(center(u), center(v))

function norm_sig_xcorr!(a, u)
    um = mean(u)
    uvar = std(u; corrected = false) * sqrt(length(u))
    a .= u .- um ./ uvar
end

function norm_sig_xcorr(u)
    a = similar(u)
    norm_sig_xcorr!(a, u)
    return a
end

xcorr_normed(u, v) = xcorr(norm_sig_xcorr(u), norm_sig_xcorr(v))

xcorr_ndx_lag(ndx::Integer, xc_len::Integer) = ndx - cld(xc_len, 2)

xcorr_lags(sig_len::Integer) = -(sig_len-1):(sig_len-1)
xcorr_lags(a::AbstractVector) = xcorr_lags(length(a))

xcorr_best_lag(xc::AbstractArray) = xcorr_ndx_lag(indmax(xc), length(xc))
