#=
sino-plot.jl
=#

# include("helper.jl")

using Plots: gui
#using MIRTjim: prompt
using Sinograms: SinoPar, SinoMoj, SinoFanArc, SinoFanFlat, sino_plot_rays

#using Unitful: mm, °
using Test: @test, @testset, @test_throws, @inferred


function _test_sino_plot(geo ; nb::Int = 60, na::Int = 40)
    d, orbit = 2, 180
#   d, orbit = 2mm, 180.0°
    @show sg = geo(; d, orbit, nb, na)
    sino_plot_rays(sg)
    true
end

@testset "sino_plot" begin
#   sg = SinoFan(Val(:ge1))
    sg_list = (SinoPar, SinoFanArc, SinoFanFlat, SinoMoj)
    for geo in sg_list
        @test _test_sino_plot(geo)
        gui()
#       prompt()
    end
end
