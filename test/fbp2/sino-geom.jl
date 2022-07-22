#=
test/fbp2/sino-geom.jl
=#

include("../helper.jl")

using Sinograms: RealU
using Sinograms: SinoGeom, SinoParallel, SinoFan
using Sinograms: SinoPar, SinoMoj, SinoFanArc, SinoFanFlat
#using Sinograms #: sino_geom, sino_geom_par, sino_geom_fan, sino_geom_moj
using Sinograms: ones, zeros, angles, rays, downsample, oversample # values
#import Sinograms: sino_geom_gamma
# todo

using Unitful: mm, °
using Test: @test, @testset, @test_throws, @inferred


function values(sg::S) where {S <: SinoGeom}
#   Iterators.map(p -> getproperty(sg, p), fieldnames(S))
    [getproperty(sg, p) for p in fieldnames(S)]
end


function _test(geo)
    sg = @NOTinferred geo(; d, orbit)
#   args = 128, 100, d, orbit, 0*orbit, 0, 2d
    args = values(sg)
    @inferred geo{Int,Float64}(args...)
    sg = @NOTinferred geo(args...)
    downsample(sg, 2)
    oversample(sg, 2)
    true
end

d = 2mm
orbit = 180.0°
# d = 2; orbit = 180

@testset "construct" begin
    sg = SinoPar(; d, orbit)
    sg = SinoMoj(; d, orbit)
    sg = SinoFanArc(; d, orbit)
    sg = SinoFanFlat(; d, orbit)
    sg = SinoFan(Val(:ge1))

    for geo in (SinoPar, SinoMoj, SinoFanArc, SinoFanFlat)
        @test _test(geo)
    end

    @test SinoFan(Val(:ge1) ; orbit=:short) isa SinoFanArc
end




# sino_geom.jl

#ig = image_geom(nx=512, fov=500).down(down)
#ig = image_geom(nx=ig.nx, dx=ig.dx, mask=ig.circ())

sg_list = (SinoPar, SinoMoj, SinoFanArc, SinoFanFlat)

#=
sg_list = (
    sino_geom(:par ; down=down, d=4),
    sino_geom(:ge1 ; down=down, orbit_start=20, dfs=0), # arc
    sino_geom(:ge1 ; down=down, orbit_start=20, dfs=Inf), # flat
    sino_geom(:moj ; down=down, d=4*sqrt(2)),
    sino_geom(:fan ; down=down, orbit=:short),
)
=#

function _test_prop(sg)
    sg.ad[2]
    sg.rfov
    @NOTinferred sg.down(2)
    @NOTinferred sg.over(2)
    sg.dim
    sg.w
    sg.ones
    sg.zeros
    sg.dr
    sg.ds
    sg.r
    sg.s
    sg.ad
    sg.ar
    sg.xds
    sg.yds

    show(isinteractive() ? stdout : devnull, sg)
    show(isinteractive() ? stdout : devnull, MIME("text/plain"), sg)
    @test ones(sg) isa AbstractMatrix{Float32}
    @test zeros(sg) isa AbstractMatrix{Float32}
    @test angles(sg) isa AbstractVector
    @test rays(sg) isa Tuple
    @test propertynames(sg) isa NTuple

    if sg isa SinoParallel
        @test isnothing(sg.gamma)
    end

    if sg isa SinoFan
        @test sg.gamma isa AbstractVector
        sg.gamma_max
        @test sg.orbit_short isa Real
        sg.dsd
        sg.dfs
        sg.dso
#       sg.dod
    end

    if sg isa SinoMoj
        sg.d_moj(0)
        sg.d_ang # angular dependent d for :moj
    end

    @test sg.shape(vec(sg.ones)) == sg.ones
    x = (1:4) * oneunit(d)
    sg.taufun(x, 2*x)
    sg.unitv()
    sg.unitv(ib=1, ia=2)
    true
end

down = 8
orbit_start = 20°
for geo in sg_list
    sg = geo( ; d, orbit, orbit_start, down)
    @test _test_prop(sg)
end


#=
    sg.grid
    sg.plot_grid(plot)
#    prompt()

nsg = length(sg_list)
pl = Array{Any}(undef, nsg)
for ii = 1:nsg

# pl[ii] = plot(); sg.plot!(plot! ; ig=ig)

#@inferred todo
@test sino_geom(:ge1 ; units=:cm) isa SinoFanArc

@test_throws String sino_geom(:badhow)
@test_throws String sino_geom(:ge1 ; dfs=-1)
@test_throws String sino_geom(:ge1 ; units=:bad)

pg = sino_geom_plot_grids(plot)
ps = sino_geom_show()

plot(pg..., ps..., layout=(2,4))
prompt()
plot(pl[1:4]...)
=#
