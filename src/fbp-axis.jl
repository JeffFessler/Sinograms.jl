#=
fbp-axis.jl
Simple interfaces to FBP for AxisArrays for user convenience
WIP
=#

#=
using .AxisArrays: AxisMatrix
using Sinograms: FBPplan, fbp!

export fbp


"""
    fbp(sino::AbstractMatrix{<:Number} ; r, ϕ, args...)

FBP reconstruction from a parallel-beam sinogram of size `[nr × nϕ]`.
Returns an image of size `[nr × nr]`.

# Options
* `r` : radial sampling; by default `-(nr-1)/2 : (nr-1)/2`
* `ϕ` : angular sampling; by default `(0:nϕ-1)/nϕ * π`
"""
function fbp(
    sino::AbstractMatrix{<:Number} ;
    r::StepRangeLen{<:Number} = LinRange(-1,1,size(sino,1)) * (size(sino,1)-1)/2,
    ϕ::StepRangeLen{<:Number} = (0:size(sino,2)-1) / size(sino,2) * π,
    kwargs...,
)
    fbp(AxisArray(sino; r, ϕ) ; kwargs...)
end


"""
    fbp(sino::AxisMatrix ; args...)

"""
function fbp(
    sino::AxisMatrix{<:Number} ;
    nx::Int = size(sino, 1),
    ny::Int = ny,
    dx::RealU = diff(axes(sino)[1][[2,1]]), # dr
    dy::RealU = dx,
    x::StepRangeLen{<:Number} = LinRange(-1,1)*(nx-1)/2/nx * dx,
    y::StepRangeLen{<:Number} = LinRange(-1,1)*(ny-1)/2/ny * dy,
)
    T = promote_type(eltype(sino), Float32)
    image = AxisMatrix(zeros(T, nx, ny), x, y)
    fbp!(image, sino)
end


function fbp!(
    image::AxisMatrix{<:Number},
    sino::AxisMatrix{<:Number},
)
    sg = sino_geom(:par, sino)
    ig = image_geom(image)
    plan = FBPplan(sg, ig)
    fbp!(image, sino, plan)
end


function sino_geom(how::Symbol, sino::AxisMatrix)
    sino_geom(how,
        nb = size(sino,1),
        na = size(sino,2),
        dr = diff(axes(sino)[1][[2,1]])
        da = diff(axes(sino)[2][[2,1]])
        orbit_start = axes(sino)[2][1]
        orbit = axes(sino)[2][1]
    )
end


function image_geom(image::AxisMatrix)
    image_geom(
        nx = size(image,1),
        ny = size(image,2),
        dx = diff(axes(image)[1][[2,1]])
        dy = diff(axes(image)[2][[2,1]])
        offset_x = axes(image)[1][1] / dx + nx/2-1
        offset_y = axes(image)[2][1] / dx + nx/2-1
    )
end
=#
