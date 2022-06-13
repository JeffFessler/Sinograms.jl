include("ct_geom.jl")
include("cylinder_proj.jl")
include("ellipsoid_proj.jl")
using MIRT
using MIRTjim
#using Plots
#using PlotlyJS
include("feldkamp.jl")
down = 4
fover = 2
cg = ct_geom(:arc ; ns = 256, nt = 240, na = 288, ds = 4, dt = 4, down = down, offset_s = 0.25, offset_t = 0, dsd = 949.075, dod = 408.075)
#println("Value of rmax ")
#println(cg.rmax)

fov = 500
nx = 64
ny = 60
nz = 50
dims = (nx, ny, nz)
delta = fov/nx.*(1, -1, 1)
mask2 = ones(nx,ny) .== ones(nx,ny)
mask2[nx,ny] = 0
mask = repeat(mask2, inner = (1,1,1), outer = (1,1,nz))
println("Image geometry")
#=
subtle difference between matlab code: matlab can handle 3D image geometry being
input into ellipse_im even though it's supposed to be 2D, julia cannot
=#
ig2 = MIRT.ImageGeom(dims = dims[1:2], deltas = delta[1:2], offsets = (0,0), mask = mask[:,:,1])
ig3 = MIRT.ImageGeom(dims = dims, deltas = delta, offsets = (0,0,0), mask = mask)

cyl = [20, 10, 10, 150, 150, Inf, 0, 0]
ell = [20 10 10 150 150 980 0 0 0.01; 80 10 10 50 50 30 0 0 0.01; -10 -40 75 40 40 40 0 0 0.01; -10 80 -20 30 30 30 0 0 0.01]
println("Ellipsoid object initialized")

fover = 2
inds = [1,2,4,5,7,8]
tmp = ellipse_im(ig2, reshape(cyl[inds], (1,6)) ; oversample = fover)
#in line 48 of feldkamp_example.m, tmp is all zeros and has two additional parameters
tmp = repeat(tmp, inner = (1,1,1), outer = (1,1,ig3.nz))
xtrue = tmp + ellipsoid_im(ig3, ell ; oversample = fover)
#jim(xtrue)
#prompt()
#plot(PlotlyJS.heatmap(z=xtrue[:,:,24]))
#Plots.heatmap(xtrue[:,:,24])

tmp = cylinder_proj(cg, reshape(cyl, (1,8)), fover)
proj = tmp + ellipsoid_proj(cg, ell, fover)
#Plots.heatmap(proj[:,:,1])

li_hat = proj

xfdk, junk = feldkamp(cg, ig3, li_hat)
errorvector = xfdk - xtrue
println(minimum(errorvector))
println(maximum(errorvector))

cg.plot3(ig3)

#=
Run the code below to look at Feldkamp results for central and offcenter slices
ix = 1:ig3.nx
iy = convert(Int,ceil(ig3.ny/2))
iz = convert(Int,ceil(ig3.nz/2))
Plots.plot(ix, xtrue[ix,iy,iz], xlim=(1,ig3.nx), ylim=(-0.003, 0.023))
Plots.plot!(ix, xfdk[ix, iy, iz])

iz = 1:ig3.nz
ix = convert(Int, 1+floor(ig3.nx/2))
iy = convert(Int, 1+floor(ig3.ny/2))
Plots.plot(iz, xtrue[ix,iy,iz], xlim=(1,ig3.nz), ylim=(-0.003, 0.023))
Plots.plot!(iz, xfdk[ix, iy, iz])
=#
