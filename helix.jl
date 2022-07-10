function rebin_helix(cg, ig, proj, varargin)

    return sino, orbits, used
end

function rebin_helix_do(proj, mask2, zslice, dz, dx, nx, 
    ny, nz, na, ns, ds, ss,
    nt, dt, ## tpoints,
    dsd, dso, dod, dfs, offset_s, wt,
    orbit, orbit_start, pitch, rmax,
    source_zs, gamma,
    type, short, chat)

    if dfs != 0 && dfs != Inf 
        return "only arc and flat done" ##Fail statement
    end

    na1 = na/(orbit*360); 
    if abs(floor(int, na1) - na1) > 1e-4
        return "bad na/orbit" ##Fail statement
    end

    orbits = zeros(nz,2)

    if short ##------------
        orbit_short_ideal = 180 + 2 * rad2deg(max(abs(gamma)))
        na1 = 2 * ceil(orbit_short_ideal / (orbit / na) / 2) 
        orbit_short = na1 * (orbit / na)
        na1_half = na1 / 2 ## because even
        orbits[:,1] .= orbit_short
    else
        orbits[:,1] .= 360
	    na1_half = ceil(na1 / 2)
    end

    sino = zeros(ns, na1, nz)
    used = zeros(nt, na)

    for iz = 1:nz ##---------
        ## ticker(mfilename, iz, nz) ????????
        zmid = zslice[iz] ## z coordinate of center of this slice
	    ia1_middle = imin[abs(zmid - source_zs)] ## "1" because julia indexing
	    ia1_list = ia1_middle - na1_half .+ collect(0:(na1-1))
        
        if ia1_list[1] < 1
            ia1_list = collect(1:na1)
        elseif ia1_list[length(ia1_list)] > na
            ia1_list = collect(1:na1) .- na1 .+ na ## Subtracting the two arrays and then adding na or vice versa
        end
        orbits[iz,2] = orbit_start + (ia1_list[1]-1) / na * orbit

        for i1=1:na1
            ia1 = ia1_list[i1]
            zdiff = zmid - source_zs[ia1]
		    if dfs == Inf
			    tt = zdiff * (dsd^2 .+ ss.^2) / (dsd * dso) ##assuming ss is array
            elseif dfs == 0
			    tt = zdiff * (dsd / dso) ./ cos.(gamma) ##assuming gamma is array
            else
			    return "bug"
		    end

            itr = tt./dt + wt ## Assuming the
		    itr = max(itr, 0)
		    itr = min(itr, nt-1)

            if type == round(ns,1)#---------
                itn = round(itr)
                it1 = 1 + itn
                it1 = max(it1, 1)
                it1 = min(it1, nt)
                tt = ((it1-1) - wt) * dt
    
                
                scale = rebin_helix_scale(dfs, dsd, ss, tt)
    
                tmp = collect(1:ns) + (it1-1)*ns + (ia1-1)*ns*nt
                view = proj(tmp)
        		## view = proj(:, it1, ia1)
    
                used(it1, ia1) = 1

            elseif type == Float64
                tt = (itr - wt) * dt
                it0 = floor(itr)
                it0 = min(it0, nt-2)
                frac = itr - it0
    
                scale = rebin_helix_scale(dfs, dsd, ss, tt);
    
                tmp = collect(1:ns) .+ it0*ns .+ (ia1-1)*ns*nt
                tmp0 = proj[tmp]
                tmp1 = proj[tmp + ns]
    			## view = (1 - frac) .* tmp0 + frac .* tmp1
                view = tmp0 + frac .* (tmp1 - tmp0)
                used[[it0+1; it0+2], ia1] .= 1
    
            else
                return "bad type"
            end

            sino[:, i1, iz] = scale .* view
        end
    end
    return sino, orbits, used
end

function rebin_helix_scale(dfs, dsd, ss, tt) 
    if dfs == Inf
        scale = sqrt.(ss.^2 .+ tt.^2 .+ dsd^2) ./ sqrt.(tt.^2 .+ dsd^2)
    elseif dfs == 0
        scale = dsd ./ sqrt.(tt.^2 .+ dsd^2)
    else
        return "bad dfs"
    end
    return scale
end