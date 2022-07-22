# test/fbp/window.jl

using Sinograms: fbp_window, Window
using Sinograms: NoWindow, Boxcar, Hamming, Hann, WindowVect
using Test: @test, @testset, @inferred


@testset "windows" begin
    shapes = (NoWindow, Boxcar, Hamming, Hann)
    for shape in shapes
        ws = @inferred Window(shape(), 0.7)
        @test ws isa Window
        win = @inferred fbp_window(ws, 8)
        @test win[1] == 1
    end

    v = 1:8
    ws = @inferred Window(WindowVect(v))
    @test ws isa Window
    win = @inferred fbp_window(ws, 8)
    @test win[1] == 5
end
