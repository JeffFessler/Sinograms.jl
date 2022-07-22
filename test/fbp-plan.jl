#=
fbp-plan.jl
=#

using Sinograms: fbp2, SinoGeom, NormalPlan
using Sinograms: SinoPar
using ImageGeoms: ImageGeom
using Test: @test, @testset, @inferred


@testset "fbp-plan" begin
    sg = SinoPar()
    ig = ImageGeom()
    plan = @inferred NormalPlan(sg, ig, Window(), ones(3,5))
    @test plan isa NormalPlan
#   plan = @inferred fbp2(sg, ig) # todo!
    plan = fbp2(sg, ig)
    @test plan isa NormalPlan
end
