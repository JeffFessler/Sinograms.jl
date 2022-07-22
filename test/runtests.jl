# runtests.jl

using Sinograms
using Test: @test, @testset, detect_ambiguities

include("helper.jl")

# todo: more

include("fbp/window.jl")

include("fbp2/sino-geom.jl")
include("fbp2/ramp.jl")

include("fbp-plan.jl")
include("fbp-par.jl")
include("zwart_powell.jl")

@testset "Sinograms" begin
    @test isempty(detect_ambiguities(Sinograms))
end
