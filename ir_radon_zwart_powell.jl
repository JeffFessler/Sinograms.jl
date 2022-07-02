## Matlab author: Seongjin Yoon, Univ. of Michigan, Jan 2015
## Translated by: Harshit Nanda


export ir_radon_zwart_powell

"""
    output = ir_radon_zwart_powell(theta, rr)
Analytic 2D Radon transform of Zwart-Powell box spline. To use this 
you call the theta matrix which represents the Zwart-Powell element
by the Box-Spline directions and the distance between ray and the point.

in
- `theta`         ray angle in radian
- `rr`            distance between the point and the ray (normalized by the pixel size)

out
- `output`        radon transform of the Zwart-Powell box spline element
"""
function ir_radon_zwart_powell(theta, rr)
    dim = size(theta)
    theta = vec(theta)    
    zeta = zeros(length(theta), 4) 

    # setting each dimension
    zeta[:,1] = cos.(theta)
    zeta[:,2] = sin.(theta)
    zeta[:,3] .= zeta[:,1] .+ zeta[:,2]
    zeta[:,4] .= zeta[:,1] .- zeta[:,2]

    # using matlab floating point value of eps(1.0)
    cond = abs.(zeta) .>= 1.1920929e-07

    N = sum(cond, dims = 2)
    N = vec(N)      
    output = BoxSp4(rr[:], zeta, cond, N) ./ factorial.(N.-1)
    output = reshape(output, dim)
    return output
end 

## BoxSp0
BoxSp0(y, N::Int) = y > 0 ? y^(N-1) : zero(y)  

## BoxSp1
function BoxSp1(y, zeta, cond, N)
    good = cond[:,1]
    output = (BoxSp0.(y .+ (0.5 .* zeta[:,1]), N) .- BoxSp0.(y .- (0.5 .* zeta[:,1]), N)) ./ zeta[:,1]
    output[.~good] .= BoxSp0.(y[.~good], N[.~good])
    return output
end
    
## BoxSp2
function BoxSp2(y, zeta, cond, N)
    good = cond[:,2]

    # Average of second column
    output = (BoxSp1(y .+ (0.5.*zeta[:,2]), zeta, cond, N) .- BoxSp1(y .- (0.5 .* zeta[:,2]), zeta, cond, N)) ./ zeta[:,2]
    output[.~good] .= BoxSp1(y[.~good], zeta[.~good,:], cond[.~good,:], N[.~good])
    return output 
end

## BoxSp3
function BoxSp3(y, zeta, cond, N)
    good = cond[:,3]

    # Average of third column
    output = (BoxSp2(y .+ 0.5 .* zeta[:,3], zeta, cond, N) .- BoxSp2(y .- 0.5 .* zeta[:,3], zeta, cond, N)) ./ zeta[:,3]
    output[.~good] .= BoxSp2(y[.~good], zeta[.~good,:], cond[.~good,:], N[.~good])
    return output
end

## BoxSp4
function BoxSp4(y, zeta, cond, N)
    good = cond[:,4]

    # Average of fourth column
    output = (BoxSp3(y .+ 0.5 .* zeta[:,4], zeta, cond, N) .- BoxSp3(y .- 0.5 .* zeta[:,4], zeta, cond, N)) ./ zeta[:,4] 
    output[.~good] .= BoxSp3(y[.~good], zeta[.~good,:], cond[.~good,:], N[.~good])
    return output
end
