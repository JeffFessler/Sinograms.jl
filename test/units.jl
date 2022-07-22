# units.jl

using Sinograms: to_radians
using Unitful: °, rad
using Test

@testset "units" begin
    ar = @inferred to_radians(180) # unitless (like degrees)
    @test ar ≈ π
    ar = @inferred to_radians(1.0 * rad) # radians
    @test ar ≈ 1.0 * rad
    ar = @inferred to_radians(180°) # degrees
    @test ar ≈ π * rad
end
