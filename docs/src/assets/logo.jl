using MIRT: sino_geom, ellipse_sino
using MIRTjim: jim
using Colors: RGB
using FileIO: save

const julia_purple  = (0.584, 0.345, 0.698)
const julia_green   = (0.22, 0.596, 0.149)
const julia_red     = (0.796, 0.235, 0.2)
colors = [ julia_purple, julia_green, julia_red, ] # from Luxor.jl

sg = sino_geom(:par, nb=200, na=181, orbit=181)

pars = [
    [70 0 10 10 0 1];
    [30 30 10 10 0 1];
    [0 -60 10 10 0 1];
]

sino = zeros(sg.dim..., size(pars,1))
logo = zeros(RGB{Float32}, sg.dim...)
for ip = 1:3
    sino[:,:,ip] = ellipse_sino(sg, pars[ip:ip,:])
    c = RGB{Float32}(colors[ip]...)
    logo[sino[:,:,ip] .> 0] .= c
end
sino = round.(Int, 255*sino/maximum(sino))
jim(sg.r, sg.ad, sino)

#save("logo.png", logo')
