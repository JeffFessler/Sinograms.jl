using FFTW
using Images
include("fbp_ramp.jl")
function fdk_fan_filter(type, n, ds, dsd, window)
    if type == "flat"
        h, nn = fbp_ramp("flat", n, ds)
    else
        h, nn = fbp_ramp("arc", n, ds, dsd)
    end
    H = real.(fft(fftshift(h)))

    if window == "ramp"
        window = ones(n)
    elseif window == "hann"
        window = ones(n) #fix this later
    else
        window = ones(n) #fix this later
    end

    H = H .* fftshift(window)
    return H
end

function fdk_filter(proj, window, dsd, dfs, ds)
    (ns, nt, na) = size(proj)
    npadh = convert(Int64, 2^ceil(log2(2*ns-1)))

    #it might be better just to use a convolution
    if isinf(dsd)
        H = fdk_fan_filter("flat", npadh, ds, [], window)
    elseif isinf(dfs)
        H = fdk_fan_filter("flat", npadh, ds, [], window)
    elseif dfs == 0
        H = fdk_fan_filter("arc", npadh, ds, dsd, window)
    end
    H = ds * H

    #something about symmetric IFFT here
    #also how to deal with N point FFT
    #proj = ifft(fft(proj, npadh, 1)) this is the matlab code
    padded = [proj; zeros(npadh-ns, nt, na)]
    proj = ifft(fft(padded, 1) .* repeat(H, inner = (1,1,1), outer = (1, nt, na)), 1)
    proj = proj[1:ns, :, :]
    return real.(proj)
end
