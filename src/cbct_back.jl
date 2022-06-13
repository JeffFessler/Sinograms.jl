using FFTW
export cbct_back
"""
cone-beam backprojector for feldkamp.jl

in
    proj  [ns nt na]     cone-beam projection views
    cg    strum          CtGeom
    ig    strum          ImageGeom

out
    img   [nx ny nz]     back projection result

Translated from cbct_back.m in MIRT
Copyright 2022-5-18 Jason Hu and Jeff Fessler, University of Michigan
"""
function cbct_back(proj, cg, ig)
    #line 143 of cbct_back.m
    ns = convert(Int64, cg.ns)
    nt = convert(Int64, cg.nt)
    na = convert(Int64, cg.na)
    ds = cg.ds
    dt = cg.dt
    offset_s = cg.offset_s
    offset_t = cg.offset_t
    offset_source = 0
    dsd = cg.dsd
    dso = cg.dso
    dfs = cg.dfs
    orbit = cg.orbit
    orbit_start = cg.orbit_start
    source_zs = zeros(na) #for some reason in matlab this array is all zeros
    if isa(ig, im_geom2)
        mask = ig.mask
    elseif isa(ig, ImageGeom)
        mask = ig.mask_or
    else
        error("Mask Not implemented yet")
    end
    nz = convert(Int64, ig.nz)
    dx = ig.dx
    dy = ig.dy
    dz = ig.dz
    offset_xyz = [ig.offset_x, ig.offset_y, ig.offset_z]
    ia_skip = 1
    scale_dang = true

    (nx, ny) = size(mask)
    betas = deg2rad.(orbit * [0:na-1;] / na .+ orbit_start)
    wx = (nx-1)/2 + offset_xyz[1];
    wy = (ny-1)/2 + offset_xyz[2];
    wz = (nz-1)/2 + offset_xyz[3];

    xc, yc = ndgrid(([0:nx-1;] .- wx) * dx, ([0:ny-1;] .- wy) * dy)
    zc = ([0:nz-1;] .- wz) * dz

    xc = xc[mask]
    yc = yc[mask]

    ws = (ns+1)/2 + offset_s
    wt = (nt+1)/2 + offset_t

    (t1, t2) = size(mask)
    img = zeros(t1, t2, nz)
    sdim = (ns+1, nt)
    proj1 = zeros(ns+1, nt)

    for iz in 1:nz
        ia_min = 1
        ia_max = na

        img2 = zeros(ns*nt-1)
        for ia in [ia_min:ia_skip:ia_max;]
            beta = betas[ia]

            #matlab code has +xc for some reaso
            x_beta = xc * cos.(beta) + yc * sin.(beta)
            y_betas = - (-xc * sin.(beta) + yc * cos.(beta)) .+ dso

            if isinf(dsd) || isinf(dso)
                mag = 0*y_betas .+ 1
            else
                mag = dsd ./ y_betas
            end

            if isinf(dfs) || isinf(dsd) || isinf(dso)
                sprime = mag .* x_beta
            elseif dfs == 0
                r_loop = x_beta .- offset_source
                sprime = dsd * atan.(r_loop, y_betas)
            end

            tprime = mag * (zc[iz] - source_zs[ia])

            bs = sprime / ds .+ ws
            bt = tprime / dt .+ wt

            bs[bs .< 1 .|| bs .> ns] .= ns+1
            bt = max.(bt, 1)
            bt = min.(bt, nt)

            is = floor.(bs)
            it = floor.(bt)

            is[is .== ns+1] .= ns
            it[it .== nt] .= nt - 1

            wr = bs - is
            wl = -wr .+ 1
            wu = bt - it
            wd = -wu .+ 1

            proj1[1:ns, :] = proj[:,:,ia]
            is = convert(Vector{Int64}, is)
            it = convert(Vector{Int64}, it)
            p1 = wl .* proj1[Base._sub2ind(sdim, is, it)] + wr .* proj1[Base._sub2ind(sdim, is .+ 1, it)]
            p2 = wl .* proj1[Base._sub2ind(sdim, is, it .+ 1)] + wr .* proj1[Base._sub2ind(sdim, is .+ 1, it .+ 1)]

            p0 = wd .* p1 + wu .* p2

            if isinf(dfs) || isinf(dsd) || isinf(dso)
                p0 = p0 .* mag .^2
            elseif dfs == 0
                p0 = p0 .* (dsd.^2) ./ (r_loop.^2 + y_betas.^2)
            end

            img2 = img2 + p0
        end

        #matlab has embed function which in this case adds a zero then reshapes
        img[:,:,iz] = reshape([img2; 0], (ns, nt))
    end

    if scale_dang
        img = (0.5 * deg2rad(abs(orbit)) / (na/ia_skip)) * img
    end
    println("Done running")
    return img
end
