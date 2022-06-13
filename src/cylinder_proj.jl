include("ir_coord_cb_arc_to_par.jl")
include("ir_coord_cb_flat_to_par.jl")
using MIRT
export cylinder_proj

#Note: apparently this doesn't work for all angles
"""
Compute set of 2d line-integral projection views of (elliptical) cylinders.
Works for these 3D geometries:
	parallel beam
	flat-detector cone-beam
	arc-detector cone-beam (3rd generator CT)

in
	cg					CtGeom
	params [ne 8]		elliptical cylinder parameters:
			[centx, centy, centz, radx, rady, zlength, angle_degrees, amplitude]

options
	oversample			oversampling factor for emulating "strips"

out
	proj	[ns nt na]	projection views

Translated from cylinder_proj.m in MIRT
Copyright 2022-5-11, Jason Hu and Jeff Fessler, University of Michigan
"""
function cylinder_proj(cg::CtGeom, params, oversample = 1)
    if isa(cg, CtFanPar)
        return ir_cylinder_proj_do(params, cg.s, cg.t, cg.ar, cg.source_zs, Inf, cg.dod, 0, oversample)
    elseif isa(cg, CtFanArc)
        return ir_cylinder_proj_do(params, cg.s, cg.t, cg.ar, cg.source_zs, cg.dso, cg.dod, 0, oversample)
    elseif isa(cg, CtFanFlat)
        return ir_cylinder_proj_do(params, cg.s, cg.t, cg.ar, cg.source_zs, cg.dso, cg.dod, Inf, oversample)
    else
        error("Type was not found!")
    end
    return 0
end

function ir_cylinder_proj_do(params, ss, tt, beta, source_zs, dso, dod, dfs, oversample)
    if size(params, 2) != 8
        println(size(params))
        error("Cylinder must have 8 parameters")
    end
    if oversample > 1
        ds = ss[2]-ss[1]
        dt = tt[2]-tt[1]
        if any(abs.(diff(ss) / ds .- 1) .> 1e-5) || any(abs.(diff(tt) / dt .- 1) .> 1e-5)
			error("Uniform spacing required for oversampling")
		end
        No = oversample
        ss = [-(No-1):2:(No-1);] / (2*No)*ds .+ ss[:]'
        tt = [-(No-1):2:(No-1);] / (2*No)*dt .+ tt[:]'
        proj = ir_cylinder_proj_work(params, ss[:], tt[:], beta, source_zs, dso, dod, dfs, 1)
        return downsample3(proj, (No, No, 1))
    else
        return ir_cylinder_proj_work(params, ss[:], tt[:], beta, source_zs, dso, dod, dfs, 1)
    end
end

function ir_cylinder_proj_work(params, ss, tt, beta, source_zs, dso, dod, dfs, oversample)
    ns = length(ss)
    nt = length(tt)
    sss, ttt = ndgrid(ss, tt)

    if isinf(dso)
        uu = sss
        vv = ttt
        azim0 = 0 * uu
        polar = 0 * uu
    elseif isinf(dfs)
        uu, vv, azim0, polar = ir_coord_cb_flat_to_par(sss, ttt, dso, dod)
    elseif dfs == 0
        uu, vv, azim0, polar = ir_coord_cb_arc_to_par(sss,ttt,dso,dod)
    else
        error("Not done")
    end

    cpolar = cos.(polar)
    spolar = sin.(polar)
    proj = zeros(ns, nt, length(beta))

    #loop over cylinders line 102 of matlab file
    for ip = 1:size(params,1)
        par = params[ip,:]
        cx = par[1]
        cy = par[2]
        cz = par[3]
        rx = par[4]
        ry = par[5]
        zh = par[6] / 2
        eang = deg2rad(par[7])
        val = par[8]

        for ib = 1:length(beta)
            az = beta[ib] .+ azim0

            cz_eff = cz - source_zs[ib]
            ushift = cx * cos.(az) + cy * sin.(az)
            vshift = (cx * sin.(az) - cy*cos.(az)) .* spolar + cz_eff * cpolar
            az = az .- eang

            p1 = (uu-ushift) .* cos.(az) + (vv-vshift) .* sin.(az) .* spolar
            p2 = (uu-ushift) .* sin.(az) - (vv-vshift) .* cos.(az) .* spolar
            p3 = (vv-vshift) .* cpolar

            e1 = -sin.(az) .* cpolar
            e2 = cos.(az) .* cpolar
            e3 = spolar

            A = e1.^2 / rx^2 + e2.^2 / ry^2
            B = p1 .* e1 / rx^2 + p2 .* e2 / ry^2
            C = p1.^2 / rx^2 + p2.^2 / ry^2 .- 1

            det = B.^2 - A .* C
            good = det .>= 0
            tmp = sqrt.(det[good])
            l0 = 0 * det
            l1 = 0 * det
            A = A[good]
            B = B[good]
            l0[good] = (-B-tmp) ./ A
            l1[good] = (-B+tmp) ./ A

            z0 = p3 + l0 .* e3
            z1 = p3 + l1.*e3

            zswap = z0 .> z1
            z0, z1 = ir_swap(z0, z1, zswap)

            zmin = max.(z0, -zh)
            zmax = min.(z1, zh)

            l_int = 0*l1
            zok = good .& (z1 .!= z0)
            tmp0 = abs.(l1-l0)
            tmp1 = tmp0 .* max.(zmax-zmin, 0)
            tmp2 = abs.(z1-z0)
            l_int[zok] = tmp1[zok] ./ tmp2[zok]

            zeq = good .& (z1 .== z0) .& (-zh .< z0) .& (z0 .< zh)
            l_int[zeq] = tmp0[zeq]

            proj[:,:,ib] = proj[:,:,ib] + val * l_int
        end
    end
    return proj
end

function ir_swap(i0, i1, swap)
    o0 = i0
    o1 = i1
    o0[swap] = i1[swap]
    o1[swap] = i0[swap]
    return o0, o1
end
