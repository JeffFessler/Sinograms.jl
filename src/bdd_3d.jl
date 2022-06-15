## Author: Sonia Minseo Kim
# Date: June 2022

## 3-D Cone Beam Branchless Distance Driven Code
using Plots: plot
using ImagePhantoms: shepp_logan, SheppLoganToft
using MIRTjim: jim
using LazyGrids: ndgrid 
using Dierckx: Spline1D

## Forward Projection
function projection(phantom::AbstractMatrix{<:T}, geo) where {T <: Number}

    DS0 = geo.DS0
    DSD = geo.DSD
    pSize = geo.pSize
    dSize = geo.dSize
    nPix = geo.nPix
    nDet = geo.nDet
    theta = geo.theta

    for angle in 1:length(theta)

        #Source position
        sourceX = -DS0*sin(theta(angle))
        sourceY = DS0*cos(theta(angle))
        sourceZ = 0

        #Detector position
        detX = (DSD-DS0)*sin(theta(angle))
        detY = -(DSD-DS0)*cos(theta(angle))
        detZ = 0

        
        



    end #angle

end #func


function overlap_kernel()


## Begin Tests 
# Geometry Definitions
geo = (DSD = 100000, DS0 = 99700, pSize = 0.5, dSize = 0.5, nPix = 256, nDet = [512, 384],
        theta = deg2rad.(0:deg:360-deg), rotCenter = [0,0,0])

 

