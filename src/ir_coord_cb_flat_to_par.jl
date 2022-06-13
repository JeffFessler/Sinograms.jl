"""
convert from cone-beam flat panel to parallel-beam coordinates.
ss and tt must have same size; dso and dsd are scalars
in
    ss,tt       s and t coordinates of sampling points
    dso,dsd     distances for the geometry
out
	uu,vv		transaxial and axial parallel-beam detector coordinates
	azim		transaxial or azimuthal angle (radians)
	polar		polar angle (radians)

Translated from ir_coord_cb_flat_to_par.m in MIRT
2022-05-12, Jason Hu, Jeff Fessler University of Michigan
"""
function ir_coord_cb_flat_to_par(ss, tt, dso, dod)
    Ds = dso
    Dd = dod
    Dc = Ds + Dd

    uu = Ds * ss ./ sqrt.(ss.^2 .+ Dc^2)
    vv = Ds * tt ./ sqrt.(ss.^2 + tt.^2 .+ Dc^2) .* Dc ./ sqrt.(ss.^2 .+ Dc^2)

    polar = -atan.(tt ./ sqrt.(ss.^2 .+ Dc^2))

    azim = atan.(ss / Dc)
    return uu, vv, azim, polar
end
