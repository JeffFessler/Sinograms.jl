function fbp_helix_stack(cg, ig, sino, orbits, varargin...) ##varargs
    function args()

    function streq(cg)
        return cg == "test"
    end 

    ##if streq(cg, 'test')
        ## run helix_example
    ##return
    ##end
    
    varargin.window = " "
    varargin.short = true
    varargin.chat = false
    arg = vararg_pair(arg, varargin)
    
    nz = ig.nz;
    ## 2d version of image geometry
    ig = image_geom("nx", ig.nx, "dx", ig.dx, "offset_x", ig.offset_x, "ny", ig.ny, "dy", ig.dy, "offset_y", ig.offset_y,"mask", ig.mask_or)
    
    if any(orbits[:,1] .!= orbits[1,1])
        return "only constant orbit done"
    end
    
    img = ig.zeros
    sg = sino_geom("fan", "dsd", cg.dsd, "dso", cg.dso,
        "ns", cg.ns, "ds", cg.ds, "offset_s", cg.offset_s,
        "orbit", orbits(1,1),
        "na", size(sino, 2)); # views in fan-beam sinogram, not helix cg.na
    
    ## trick: Parker weights do not depend on orbit_start

    ##if varargin.short == true
        ## [parker_wt_scale180] = fbp_fan_short_wt(sg);
    ##else
        parker_wt = 1;
        scale180 = 1;
    end
    
    for iz = 1:nz   #each slice
        ticker(mfilename, iz, nz)
    
        sg.orbit_start = orbits[iz,2];
    
        ## apply Parker weighting and scaling
        sino[:,:,iz] .= scale180 * parker_wt .* sino[:,:,iz];
        geom = fbp2[sg, ig];
        img[:,:,iz] = fbp2(sino(:,:,iz), geom, "window", arg.window);
    end

    return img, sino
end
end
    