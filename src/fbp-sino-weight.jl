# fbp-sino-weight.jl

export fbp_sino_weight

using Sinograms: SinoFan, sino_geom_gamma


"""
    fbp_sino_weight(sg::SinoFan)
Return 1D sinogram weighting for first step of 2D fan-beam FBP.
"""
function fbp_sino_weight(sg::SinoFan)
    gam = sino_geom_gamma(sg)
    return @. abs(sg.dso * cos(gam) - sg.source_offset * sin(gam)) / sg.dsd
end
