"""
Hanning window with 0 only in the first element
Useful with FFT
Equivalent to matlab hann(M, 'periodic')
"""
function ir_hann_periodic(M)
    w = 0.5 * (1 .- cos.(2*pi*(0:M-1)/M))
    return w
end
