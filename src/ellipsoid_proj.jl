include("ir_coord_cb_arc_to_par.jl")
include("ir_coord_cb_flat_to_par.jl")
using MIRT
export ellipsoid_proj

"""
Compute set of 2d line-integral projection views of ellipsoid(s).
Works for both parallel-beam and cone-beam geometry.

in
	cg			ct_geom()
	params [ne 9]		ellipsoid parameters:
			[x_center y_center z_center  x_radius y_radius z_radius
				xy_angle_degrees z_angle_degrees  amplitude]
options
	oversample		over-sampling factor for emulating "strips"
				(to account for finite detector size)

out
	proj	[ns nt na]	projection views

Translated from ellipsoid_proj.m in MIRT
Copyright 2022-05-18, Jason Hu and Jeff Fessler, University of Michigan
"""
function ellipsoid_proj(cg::CtGeom, params, oversample = 1)
	if isa(cg, CtFanPar)
		return ellipsoid_proj_do(params, cg.s, cg.t, cg.ar, cg.source_zs, Inf, cg.dod, 0, oversample)
	elseif isa(cg, CtFanArc)
		return ellipsoid_proj_do(params, cg.s, cg.t, cg.ar, cg.source_zs, cg.dso, cg.dod, 0, oversample)
	elseif isa(cg, CtFanFlat)
		return ellipsoid_proj_do(params, cg.s, cg.t, cg.ar, cg.source_zs, cg.dso, cg.dod, Inf, oversample)
	else
		error("Type not acceptable")
	end
	return 0
end

function ellipsoid_proj_do(params, ss, tt, beta, source_zs, dso, dod, dfs, oversample)
    if size(params, 2) != 9
        error("Need 9 parameters per ellipsoid")
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
        proj = ellipsoid_proj_do(params, ss[:], tt[:], beta, source_zs, dso, dod, dfs, 1)
        return downsample3(proj, (No, No, 1))
    end

    ns = length(ss)
    nt = length(tt)
    sss, ttt = ndgrid(ss, tt)

    if isinf(dso)
        uu = sss
        vv = ttt
        azim0 = 0
        polar = 0
    elseif isinf(dfs)
        uu, vv, azim0, polar = ir_coord_cb_flat_to_par(sss, ttt, dso, dod)
    elseif dfs == 0
        uu, vv, azim0, polar = ir_coord_cb_arc_to_par(sss, ttt, dso, dod)
    else
        error("Not done")
    end

    cpolar = cos.(polar)
    spolar = sin.(polar)
    proj = zeros(ns, nt, length(beta))

    for ip = 1:size(params,1)
        par = params[ip,:]

        cx = par[1]
        cy = par[2]
        cz = par[3]
        rx = par[4]
        ry = par[5]
        rz = par[6]
        xang = deg2rad(par[7])
        zang = deg2rad(par[8])

        if zang != 0
            error("Z angle not done")
        end
        val = par[9]

        for ib = 1:length(beta)
            az = beta[ib] .+ azim0

            cz_eff = cz - source_zs[ib];
			ushift = cx * cos.(az) + cy * sin.(az);
			vshift = (cx * sin.(az) - cy * cos.(az)) .* spolar + cz_eff * cpolar;

			az = az .- xang;
			p1 = (uu-ushift) .* cos.(az) + (vv-vshift) .* sin.(az) .* spolar;
			p2 = (uu-ushift) .* sin.(az) - (vv-vshift) .* cos.(az) .* spolar;
			p3 = (vv-vshift) .* cpolar;

			e1 = -sin.(az) .* cpolar;
			e2 = cos.(az) .* cpolar;
			e3 = spolar;

			A = e1.^2 / rx^2 + e2.^2 / ry^2 + e3.^2 / rz^2;
			B = p1 .* e1 / rx^2 + p2 .* e2 / ry^2 + p3 .* e3 / rz^2;
			C = p1.^2 / rx^2 + p2.^2 / ry^2 + p3.^2 / rz^2 .- 1;

			proj[:,:,ib] = proj[:,:,ib] + 2 * val * real.(sqrt.(Complex.(B.^2 - A.*C))) ./ A;
        end
    end
	return proj
end

# function ellipsoid_proj_test
# 	ell = [20 0 -40 200 100 50 90 0 10; 0 50 100 80 80 20 0 0 10]
#
# end
