function ramp_arc(n, ds, dsd)
    if n/2 * ds / dsd > 0.9 * pi/2
        print("Angle is too large")
        return 0
    end
    nn = [-(n/2):(n/2-1);]
    h = 0*nn
    h[nn.==0] .= 1 / (4 * ds^2)
    odd = mod.(nn,2) .== 1
    h[odd] .= -1 ./ (pi * dsd * sin.(nn[odd] * ds / dsd)).^2
    return h, nn
end

function ramp_flat(n, ds)
    nn = [-(n/2):(n/2-1);]
    h = 0*nn
    h[convert(Int64, n/2+1)] = 1 / 4
    odd = mod.(nn,2) .== 1
    h[odd] .= -1 ./ (pi * nn[odd]).^2
    h = h / ds^2
    return h, nn
end

function fbp_ramp(type, n, ds, dsd=0)
    if mod(n,1) == 1
        print("Error, n must be even")
        return 0
    end

    if type == "arc"
        h, nn = ramp_arc(n, ds, dsd)
    elseif type == "flat"
        h, nn = ramp_flat(n, ds)
    else
        print("bad fan type")
    end
end

nb = 256
ds = 1
dsd = 100
h1, nn = fbp_ramp("arc", nb, ds, dsd)
h2, nn = fbp_ramp("flat", nb, ds)
