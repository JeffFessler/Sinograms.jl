"""
convert from cone-beam arc (3rd gen CT) to parallel-beam coordinates.
ss and tt must have same size; dso and dsd are scalars
in
    ss,tt       s and t coordinates of sampling points
    dso,dsd     distances for the geometry
out
	uu,vv		transaxial and axial parallel-beam detector coordinates
	azim		transaxial or azimuthal angle (radians)
	polar		polar angle (radians)

Translated from ir_coord_cb_arc_to_par.m in MIRT
2022-05-12, Jason Hu, Jeff Fessler, University of Michigan
"""
function ir_coord_cb_arc_to_par(ss, tt, dso, dod)
    Ds = dso
    Dd = dod
    Dc = Ds + Dd
    dfs = 0

    pos_src = [0, dso, 0]
    Rf = dfs + Dc
    pos_det = [tt;;; tt;;; tt]
    pos_det[:,:,1] = Rf * sin.(ss / Rf)
    pos_det[:,:,2] = Ds + dfs .- Rf * cos.(ss / Rf)
    #pos_det[:,:,3] = tt

    ee1 = pos_det[:,:,1] .- pos_src[1]
    ee2 = pos_det[:,:,2] .- pos_src[2]
    ee3 = pos_det[:,:,3] .- pos_src[3]
    enorm = sqrt.(ee1.^2 + ee2.^2 + ee3.^2)
    ee1 = ee1 ./ enorm
    ee2 = ee2 ./ enorm
    ee3 = ee3 ./ enorm

    polar = -asin.(ee3)
    azim = -atan.(ee1 ./ ee2)
    uu = cos.(azim) * pos_src[1] + sin.(azim) * pos_src[2]
    vv = (sin.(azim) * pos_src[1] - cos.(azim) * pos_src[2]) .* sin.(polar) + pos_src[3] * cos.(polar)
    return uu, vv, azim, polar
end
