# test/fbp2/ramp.jl


using Sinograms: SinoFanFlat, SinoFanArc, SinoPar, fbp_ramp
import Sinograms as SG
using Unitful: mm
using Test: @test, @testset, @test_throws, @inferred

include("../helper.jl")

@testset "ramp-unit" begin
    ds, dsd = 2.0f0, 99
    ds, dsd = 2.0f0mm, 99mm
    for n in 0:2
        f = @inferred SG._ramp_flat(n, ds)
        @test oneunit(f) == oneunit(1/ds^2)
        a = @inferred SG._ramp_arc(n, ds, dsd)
        @test oneunit(a) == oneunit(1/ds^2)
    end


    n = 0:2
    fun = n -> SG._ramp_flat.(n, ds)
    fun(n)
    f = @NOTinferred fun(n)
    fun = n -> SG._ramp_arc.(n, ds, dsd)
    a = @NOTinferred fun(n)

    N = 6
    ha, na = @NOTinferred SG.ramp_arc(N, ds, dsd)
    hf, nf = @NOTinferred SG.ramp_flat(N, ds)
end


@testset "fbp-ramp" begin
    N = 10
    for shape in (SinoPar, SinoFanArc, SinoFanFlat)
        sg = shape()
        h, n = @NOTinferred fbp_ramp(sg, N)
        @test h isa Vector
        @test length(h) == N
        @test length(n) == N
    end
end
