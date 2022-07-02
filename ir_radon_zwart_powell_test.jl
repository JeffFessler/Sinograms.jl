using LazyGrids: ndgrid
using Test: @test, @testset, @test_throws

@testset "Image reconstruction" begin

    ## Begin test
    theta = LinRange(0,pi,181)
    r = LinRange(-1,1,101) * 2

    # Allocating arrays using ndgrid
    (tt, rr) = ndgrid(theta, r)
    sino = ir_radon_zwart_powell(tt, rr)

    # testing if a matrix is produced
    @test sino isa Matrix
end