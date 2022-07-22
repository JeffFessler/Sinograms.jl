"""
    Sinograms module
"""
module Sinograms

using Requires: @require

const RealU = Number # Real or Unitful

include("fbp/window.jl")

include("fbp2/sino-geom.jl")
include("fbp2/ramp.jl")
include("fbp2/sino-filter.jl")

include("fbp2.jl")
include("fbp2_back.jl")
include("fbp2_back_fan.jl")
include("fbp-sino-weight.jl")

include("fbp-par.jl")
include("zwart_powell.jl")


"""
    to_radians(angle::Real)
If no Unitful package loaded, assume `angle` is in degrees and convert to radians.
"""
to_radians(angle::Real) = deg2rad(angle)

# support Plots iff user has loaded that package
function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("sino-plot.jl")
    @require Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d" include("units.jl")
end


end # module
