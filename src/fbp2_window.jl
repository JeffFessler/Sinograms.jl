function my_boxcar(n, cutoff)
    w = round(cutoff * n)
    ii = [0:n-1;] .- n/2
    window = abs.(ii) .< w/2
    return window
end

function my_hann(n, cutoff)
    w = round(cutoff * n)
    ii = [0:n-1;] .- n/2
    window = 0.5 * (1 .+ cos.(2*pi*ii/w)) .* (abs.(ii) .< w/2)
    return window
end

function my_hamming(n, cutoff)
    w = round(cutoff * n)
    ii = [0:n-1;] .- n/2
    window = (0.54 .+ 0.46 * cos.(2*pi*ii/w)) .* (abs.(ii) .< w/2)
    return window
end

function sscanf(thing::String)
    arr = split(thing, " ")
    intarr = parse.(Float32, arr)
    return intarr
    error("not done yet")
    return 0
end

function fbp2_window(n::Int, window)
    if typeof(window) == String
        if window == "" || window == "boxcar" || window == "ramp"
            window = ones(n,1)
        elseif window[1:7] == "boxcar,"
            cut = sscanf(window[8:end])
            window = my_boxcar(n, cut)
        elseif window[1:8] == "hamming,"
            cut = sscanf(window[9:end])
            window = my_hamming(n, cut)
        elseif window == "hann"
            window = my_hann(n, 1.0)
        elseif window == "hann50"
            window = my_hann(n, 0.5)
        elseif window == "hann75"
            window = my_hann(n, 0.75)
        elseif window == "hann80"
            window = my_hann(n, 0.8)
        else
            error("unknown window")
        end
    elseif size(window) != n
        error("bad window length")
    end

    window = fftshift(window)
    return window
end
