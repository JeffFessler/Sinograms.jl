#=
sino-geom.jl
sinogram geometry definitions for 2D tomographic image reconstruction
2019-07-01, Jeff Fessler, University of Michigan
2022-01-22, copied from MIRT.jl here and updated to support Unitful values
=#

using Sinograms: RealU

export SinoGeom, SinoParallel, SinoFan
export sino_geom, sino_geom_par, sino_geom_fan, sino_geom_moj
export SinoPar, SinoMoj, SinoFanArc, SinoFanFlat
export sino_geom_help


"""
    SinoGeom

Abstract type for representing sinogram ray geometries

Derived values (available by `getproperty`), i.e., `sg.?`

* `.dim`       dimensions: `(nb,na)`
* `.ds|dr`     radial sample spacing (`NaN` for `:moj`)
* `.s`         [nb] s sample locations
* `.w`         `(nb-1)/2 + offset` ('middle' sample position)
* `.ad`        source angles [degrees]
* `.ar`        source angles [radians]
* `.ones`      `ones(Float32, nb,na)`
* `.zeros`     `zeros(Float32, nb,na)`
* `.rfov`      radial fov
* `.xds`       [nb] center of detector elements (beta=0)
* `.yds`       [nb] ''
* `.grid`      (rg, phigrid) [nb na] parallel-beam coordinates # todo: cut
* `.rays`      `(rg, phigrid)` [nb na] parallel-beam coordinates
* `.plot_grid(scatter)`    plot `sg.grid` using `Plots.scatter`

For mojette:

* `.d_ang`       [na]

For fan beam:

* `.gamma`       [nb] gamma sample values [radians]
* `.gamma_max`   half of fan angle [radians]
* `.dso`         dsd - dod, Inf for parallel beam

Methods

* `.down(down)`     reduce sampling by integer factor
* `.shape(sino)`    reshape sinograms into array [nb na :]
* `.unitv(;ib,ia)`  unit 'vector' with single nonzero element
* `.taufun(x,y)`    projected s/ds for each (x,y) pair [numel(x) na]
* `.plot!(plot!;ig)`plot system geometry (mostly for SinoFan)
"""
abstract type SinoGeom end

abstract type SinoParallel <: SinoGeom end
abstract type SinoFan <: SinoGeom end


# promote_type version that checks for compatible units
function _promoter(xs::RealU...)
    all(==(oneunit(xs[1])), oneunit.(xs)) || throw("incompatible units")
    T = promote_type(eltype.(xs)...)
    T = promote_type(T, eltype(1.0f0 * oneunit(T)))
end


"""
    SinoPar : 2D parallel-beam sinogram geometry

```jldoctest
julia> SinoPar()
SinoPar{Float32, Float32} :
 nb::Int64 128
 na::Int64 200
 d::Float32 1.0
 orbit::Float32 180.0
 orbit_start::Float32 0.0
 offset::Float32 0.0
 strip_width::Float32 1.0
```
"""
struct SinoPar{Td,To} <: SinoParallel
    nb::Int           # # of "radial" samples, aka nr
    na::Int           # # of angular samples, aka nϕ
    d::Td             # dr, "radial" sample spacing
    orbit::To         # [degrees] by default, or Unitful
    orbit_start::To   # [degrees] by default, or Unitful
    offset::Float32   # sample offset, cf offset_r or offset_s [unitless]
    strip_width::Td   # same units as `d`

#=
    function SinoPar(
        nb::Int,
        na::Int,
        d::RealU,
        orbit::RealU,
        orbit_start::RealU,
        offset::Real,
        strip_width::RealU,
    )
        To = _promoter(orbit, orbit_start)
        Td = _promoter(d, strip_width)
        return new{Td,To}(nb, na,
            Td(d), To(orbit), To(orbit_start), Float32(offset), Td(strip_width))
    end
=#

    # constructor with positional arguments requires appropriate matched types
    function SinoPar(
        nb::Int,
        na::Int,
        d::Td,
        orbit::To,
        orbit_start::To,
        offset::Real,
        strip_width::Td,
    ) where {Td <: RealU, To <: RealU}
        return new{Td,To}(nb, na, d, orbit, orbit_start, Float32(offset), strip_width)
    end

end

SinoPar{Td,To}(args...) where {Td,To} = SinoPar(args...)

"""
    SinoPar( ; nb, na, d=1, orbit=180, orbit_start=0, offset=0, strip_width=d, down=0)

Constructor with named keywords.
* `orbit` and `orbit_start` must both be unitless (degrees) or have same units.
* `d` and `strip_width` must both be unitless or have the same units.

* `nb::Int = 128` # of radial bins
* `na::Int = 2 * floor(Int, nb * π/2 / 2)` # of projection angles
"""
function SinoPar( ;
    nb::Int = 128,
    na::Int = 2 * floor(Int, nb * π/2 / 2),
    d::RealU = 1,
    orbit::RealU = 180,
    orbit_start::RealU = zero(eltype(orbit)),
    offset::Real = 0,
    strip_width::RealU = d,
    down::Int = 1,
)
    To = _promoter(orbit, orbit_start)
    Td = _promoter(d, strip_width)
    sg = SinoPar(nb, na, #d, orbit, orbit_start, offset, strip_width
            Td(d), To(orbit), To(orbit_start), offset, Td(strip_width))
    return down == 1 ? sg : downsample(sg, down)
end



"""
    SinoMoj : 2D Mojette sinogram geometry

```jldoctest
julia> SinoMoj()
SinoMoj{Float32, Float32} :
 nb::Int64 128
 na::Int64 200
 d::Float32 1.0
 orbit::Float32 180.0
 orbit_start::Float32 0.0
 offset::Float32 0.0
 strip_width::Float32 1.0
```
"""
struct SinoMoj{Td,To} <: SinoParallel
    nb::Int         # # of "radial" samples, aka ns
    na::Int         # # of angular samples, aka nϕ
    d::Td           # dx, pixels must be square
    orbit::To       # [degrees]
    orbit_start::To # [degrees]
    offset::Float32 # sample offset, cf offset_r or offset_s [unitless]
    strip_width::Td #

    # constructor with positional arguments requires appropriate matched types
    function SinoMoj(
        nb::Int,
        na::Int,
        d::Td,
        orbit::To,
        orbit_start::To,
        offset::Real,
        strip_width::Td,
    ) where {Td <: RealU, To <: RealU}
        return new{Td,To}(nb, na, d, orbit, orbit_start, Float32(offset), strip_width)
    end

end

SinoMoj{Td,To}(args...) where {Td,To} = SinoMoj(args...)

"""
    SinoMoj( ; nb, na, d=1, orbit=180, orbit_start=0, offset=0, strip_width=d, down=0)

Constructor with named keywords.
* `orbit` and `orbit_start` must both be unitless (degrees) or have same units.
* `d` and `strip_width` must both be unitless or have the same units.

* `nb::Int = 128` # of radial bins
* `na::Int = 2 * floor(Int, nb * π/2 / 2)` # of projection angles
"""
function SinoMoj( ;
    nb::Int = 128,
    na::Int = 2 * floor(Int, nb * π/2 / 2),
    d::RealU = 1,
    orbit::RealU = 180,
    orbit_start::RealU = zero(eltype(orbit)),
    offset::Real = 0,
    strip_width::RealU = d,
    down::Int = 1,
)
    To = _promoter(orbit, orbit_start)
    Td = _promoter(d, strip_width)
    sg = SinoMoj(nb, na, #d, orbit, orbit_start, offset, strip_width
            Td(d), To(orbit), To(orbit_start), offset, Td(strip_width))
    return down == 1 ? sg : downsample(sg, down)
end


"""
   SinoFanArc : 2D Fan-beam sinogram geometry for arc detector

```jldoctest
julia> SinoFanArc()
SinoFanArc{Float32, Float32} :
 nb::Int64 128
 na::Int64 200
 d::Float32 1.0
 orbit::Float32 360.0
 orbit_start::Float32 0.0
 offset::Float32 0.0
 strip_width::Float32 1.0
 source_offset::Float32 0.0
 dsd::Float32 512.0
 dod::Float32 128.0
 dfs::Float32 0.0
```
"""
struct SinoFanArc{Td,To} <: SinoFan
    nb::Int           # # of "radial" samples, aka ns
    na::Int           # # of angular samples, aka nβ
    d::Td             # ds detector sample spacing
    orbit::To         # [degrees]
    orbit_start::To   # [degrees]
    offset::Float32   # sample offset, cf offset_r or offset_s [unitless]
    strip_width::Td   #
    source_offset::Td # same units as d, etc., e.g., [mm] (use with caution!)
    dsd::Td           # dis_src_det, Inf for parallel beam
    dod::Td           # dis_iso_det
#   dso::Td           # dis_src_iso = dsd-dod, Inf for parallel beam
    dfs::Td           # distance from focal spot to source

    function SinoFanArc( ;
        nb::Int = 128,
        na::Int = 2 * floor(Int, nb * π/2 / 2),
        d::RealU = 1,
        orbit::RealU = 360,
        orbit_start::RealU = zero(eltype(orbit)),
        offset::Real = 0,
        strip_width::RealU = d,
        source_offset::RealU = zero(eltype(d)),
        dsd::RealU = 4 * nb * d,
        dod::RealU = nb * d,
        dfs::RealU = zero(eltype(d)),
    )
        To = _promoter(orbit, orbit_start)
        Td = _promoter(d, strip_width, source_offset, dsd, dod, dfs)
        return new{Td,To}(nb, na,
            Td(d), To(orbit), To(orbit_start), Float32(offset), Td(strip_width),
            Td(source_offset), Td(dsd), Td(dod), Td(dfs))
    end
end


"""
    SinoFanFlat : 2D Fan-beam sinogram geometry for flat detector

```jldoctest
julia> SinoFanFlat()
SinoFanFlat{Float64, Float32} :
 nb::Int64 128
 na::Int64 200
 d::Float64 1.0
 orbit::Float32 360.0
 orbit_start::Float32 0.0
 offset::Float32 0.0
 strip_width::Float64 1.0
 source_offset::Float64 0.0
 dsd::Float64 512.0
 dod::Float64 128.0
 dfs::Float64 Inf
```
"""
struct SinoFanFlat{Td,To} <: SinoFan
    nb::Int           # # of "radial" samples, aka ns
    na::Int           # # of angular samples, aka nβ
    d::Td             # ds detector sample spacing
    orbit::To         # [degrees]
    orbit_start::To   # [degrees]
    offset::Float32   # sample offset, cf offset_r or offset_s [unitless]
    strip_width::Td   #
    source_offset::Td # same units as d, etc., e.g., [mm] (use with caution!)
    dsd::Td           # dis_src_det, Inf for parallel beam
    dod::Td           # dis_iso_det
#   dso::Td           # dis_src_iso = dsd-dod, Inf for parallel beam
    dfs::Td           # distance from focal spot to source

    function SinoFanFlat( ;
        nb::Int = 128,
        na::Int = 2 * floor(Int, nb * π/2 / 2),
        d::RealU = 1,
        orbit::RealU = 360,
        orbit_start::RealU = zero(eltype(orbit)),
        offset::Real = 0,
        strip_width::RealU = d,
        source_offset::RealU = zero(eltype(d)),
        dsd::RealU = 4 * nb * d,
        dod::RealU = nb * d,
        dfs::RealU = Inf * oneunit(d),
    )
        To = _promoter(orbit, orbit_start)
        Td = _promoter(d, strip_width, source_offset, dsd, dod, dfs)
        return new{Td,To}(nb, na,
            Td(d), To(orbit), To(orbit_start), Float32(offset), Td(strip_width),
            Td(source_offset), Td(dsd), Td(dod), Td(dfs))
    end
end


#=
"""
    function sg = sino_geom(...)

Constructor for `SinoGeom`

Create the "sinogram geometry" structure that describes the sampling
characteristics of a given sinogram for a 2D parallel or fan-beam system.
Using this structure facilitates "object oriented" code.
(Use `ct_geom()` instead for 3D axial or helical cone-beam CT.)

in
- `how::Symbol`    `:fan` (fan-beam) | `:par` (parallel-beam) | `:moj` (mojette)

options for all geometries (including parallel-beam):
- `units::Symbol`    e.g. `:cm` or `:mm`; default: :none
- `orbit_start`        default: 0
- `orbit`            [degrees] default: `180` for parallel / mojette
and `360` for fan
   * can be `:short` for fan-beam short scan
- `down::Int`        down-sampling factor, for testing

- `nb`                # radial samples cf `nr` (i.e., `ns` for `:fan`)
- `na`                # angular samples (cf `nbeta` for `:fan`)
- `d`                radial sample spacing; cf `dr` or `ds`; default 1
  * for mojette this is actually `dx`
- `offset`            cf `offset_r` `channel_offset` unitless; default 0
   * (relative to centerline between two central channels).
   * Use 0.25 or 1.25 for "quarter-detector offset"
- `strip_width`        detector width; default: `d`

options for fan-beam
- `source_offset`        same units as d; use with caution! default 0
fan beam distances:
- `dsd` cf `dis_src_det` default: `Inf` (parallel beam)
- `dod` cf `dis_iso_det` default: `0`
- `dfs` cf `dis_foc_src` default: `0` (3rd generation CT arc),
   * use `Inf` for flat detector

out
- `sg::SinoGeom`    initialized structure

See also
- `sino_geom_help()` help on methods
- `sino_geom_plot_grids()` show sampling

Jeff Fessler, University of Michigan
"""
function sino_geom(how::Symbol ; kwarg...)
    if how === :par
        sg = sino_geom_par( ; kwarg...)
    elseif how === :fan
        sg = sino_geom_fan( ; kwarg...)
    elseif how === :moj
        sg = sino_geom_moj( ; kwarg...)
    elseif how === :ge1
        sg = sino_geom_ge1( ; kwarg...)
#=
    elseif how === :hd1
        sg = sino_geom_hd1( ; kwarg...)
    elseif how === :revo1fan
        tmp = ir_fan_geom_revo1(type)
        sg = sino_geom(:fan, tmp{:}, varargin{:})
=#
    else
        throw("unknown sino type $how")
    end

    return sg
end

SinoFanMaker(dfs::Real) =
    (dfs == 0) ? SinoFanArc :
    isinf(dfs) ? SinoFanFlat :
    throw("dfs $dfs") # must be 0 or Inf
=#


# common to all
function _downsample(sg::SinoGeom, down::Int)
    nb = 2 * max(sg.nb ÷ 2down, 1) # keep it even
    na = max(sg.na ÷ down, 1)
    return (# sg.units,
        nb, na, sg.d * down, sg.orbit, sg.orbit_start, sg.offset,
        sg.strip_width * down)
end


"""
    sg = downsample(sg, down)
down-sample (for testing with small problems)
"""
function downsample(sg::T, down::Int) where {T <: SinoParallel}
    return (down == 1) ? sg : T(_downsample(sg, down)...)
end

function downsample(sg::SinoFan, down::Int)
    down == 1 && return sg
    return SinoFanMaker(sg.dfs)(_downsample(sg, down)...,
        sg.source_offset, sg.dsd, sg.dod, sg.dfs)
end


#=
# common to all
function _sino_geom_over(sg::SinoGeom, over::Int)
    return (sg.units,
        sg.nb * over, sg.na, sg.d / over,
        sg.orbit, sg.orbit_start, sg.offset * over,
        sg.strip_width / over)
end

"""
    sg = sino_geom_over(sg, over::Int)
over-sample in "radial" dimension
For Mojette sampling, it means that `d = dx/over`.
"""
function sino_geom_over(sg::T, over::Int) where {T <: SinoParallel}
    return (over == 1) ? sg : T(_sino_geom_over(sg, over)...)
end
function sino_geom_over(sg::SinoFan, over::Int)
    return (over == 1) ? sg :
        SinoFanMaker(sg.dfs)(_sino_geom_over(sg, over)...,
            sg.source_offset, sg.dsd, sg.dod, sg.dfs)
end


"""
    sg = sino_geom_fan()
"""
function sino_geom_fan( ;
    units::Symbol = :none,
    nb::Int = 128,
    na::Int = 2 * floor(Int, nb * π/2 / 2),
    d::Real = 1,
    orbit::Union{Symbol,Real} = 360, # [degrees]
    orbit_start::Real = 0,
    strip_width::Real = d,
    offset::Real = 0,
    source_offset::Real = 0,
    dsd::Real = 4*nb*d,    # dis_src_det
#    dso::Real = [],        # dis_src_iso
    dod::Real = nb*d,    # dis_iso_det
    dfs::Real = 0,        # dis_foc_src (3rd gen CT)
    down::Int = 1,
)

    maker = (dfs == 0) ? SinoFanArc :
            isinf(dfs) ? SinoFanFlat :
            throw("dfs $dfs") # must be 0 or Inf

    if orbit === :short # trick
        sg_tmp = maker(units,
            nb, na, d, 0, orbit_start, offset, strip_width,
            source_offset, dsd, dod, dfs)
        orbit = sg_tmp.orbit_short
    end
    isa(orbit, Symbol) && throw("orbit :orbit")

    sg = maker(units,
        nb, na, d, orbit, orbit_start, offset, strip_width,
        source_offset, dsd, dod, dfs)

    return downsample(sg, down)
end


"""
    sg = sino_geom_par( ... )
"""
function sino_geom_par( ;
    units::Symbol = :none,
    nb::Int = 128,
    na::Int = 2 * floor(Int, nb * π/2 / 2),
    down::Int = 1,
    d::Real = 1,
    orbit::Real = 180, # [degrees]
    orbit_start::Real = 0,
    strip_width::Real = d,
    offset::Real = 0,
)
    sg = SinoPar(units,
        nb, na, d, orbit, orbit_start, offset, strip_width)
    return downsample(sg, down)
end


"""
    sg = sino_geom_moj( ... )
"""
function sino_geom_moj( ;
    units::Symbol = :none,
    nb::Int = 128,
    na::Int = 2 * floor(Int, nb * π/2 / 2),
    down::Int = 1,
    d::Real = 1, # means dx for :moj
    orbit::Real = 180, # [degrees]
    orbit_start::Real = 0,
    strip_width::Real = d, # ignored ?
    offset::Real = 0,
)
    sg = SinoMoj(units,
        nb, na, d, orbit, orbit_start, offset, strip_width)
    return downsample(sg, down)
end


"gamma for general finite dfs (rarely used)"
function sino_geom_gamma_dfs(sg::SinoFan)
    dis_foc_det = sg.dfs + sg.dsd
    alf = sg.s / dis_foc_det
    atan.(dis_foc_det * sin.(alf), dis_foc_det * cos.(alf) .- sg.dfs)
    # equivalent to s/dsd when dfs=0
end


"""
    sino_geom_gamma()
gamma sample values for :fan
"""
function sino_geom_gamma(sg::SinoFan)
    return    sg.dfs == 0 ? sg.s / sg.dsd : # 3rd gen: equiangular
            isinf(sg.dfs) ? atan.(sg.s / sg.dsd) : # flat
            sino_geom_gamma_dfs(sg) # general
end


"""
    sino_geom_rfov()
radial FOV
"""
sino_geom_rfov(sg::SinoPar) = maximum(abs.(sg.r))
sino_geom_rfov(sg::SinoFan) = sg.dso * sin(sg.gamma_max)
sino_geom_rfov(sg::SinoMoj) = sg.nb/2 * minimum(sg.d_ang) # (ignores offset)


function _sino_geom_taufun(sg::SinoParallel, x, y)
    ar = sg.ar' # row vector, for outer-product
    return (x * cos.(ar) + y * sin.(ar)) / sg.dr # tau
end

function _sino_geom_taufun(sg::SinoFan, x, y)
    b = sg.ar' # row vector, for outer-product
    xb = x * cos.(b) + y * sin.(b)
    yb = -x * sin.(b) + y * cos.(b)
    tangam = (xb .- sg.source_offset) ./ (sg.dso .- yb) # e,tomo,fan,L,gam
    if sg.dfs == 0 # arc
        tau = sg.dsd / sg.ds * atan.(tangam)
    elseif isinf(sg.dfs) # flat
        tau = sg.dsd / sg.ds * tangam
#    else
#        throw("bad dfs $(sg.dfs)")
    end
    return tau
end


"""
    sino_geom_taufun()
projected `s/ds`, useful for footprint center and support
"""
function sino_geom_taufun(sg::SinoGeom, x, y)
    size(x) != size(y) && throw("bad x,y size")
    return _sino_geom_taufun(sg, vec(x), vec(y))
end


"""
    sino_geom_xds()
center positions of detectors (for beta = 0)
"""
sino_geom_xds(sg::SinoPar) = sg.s
sino_geom_xds(sg::SinoMoj) = sg.s # todo: really should be angle dependent
sino_geom_xds(sg::SinoFanArc) = sg.dsd * sin.( sg.gamma) .+ sg.source_offset
sino_geom_xds(sg::SinoFanFlat) = sg.s .+ sg.source_offset


"""
    sino_geom_yds()
center positions of detectors (for beta = 0)
"""
sino_geom_yds(sg::SinoPar) = zeros(Float32, sg.nb)
sino_geom_yds(sg::SinoMoj) = zeros(Float32, sg.nb)
sino_geom_yds(sg::SinoFanArc) = sg.dso .- sg.dsd * cos.(sg.gamma)
sino_geom_yds(sg::SinoFanFlat) = fill(-sg.dod, sg.nb)


"""
    sino_geom_unitv()
sinogram with a single ray
"""
function sino_geom_unitv(
    sg::SinoGeom ;
    ib::Int = sg.nb ÷ 2 + 1,
    ia::Int = sg.na ÷ 2 + 1,
)
    out = sg.zeros
    out[ib,ia] = 1
    return out
end


"""
    (rg, ϕg) = sino_geom_grid(sg::SinoGeom)

Return grids `rg` and `ϕg` (in radians) of size `[nb na]`
of equivalent *parallel-beam* `(r,ϕ)` (radial, angular) sampling positions,
for any sinogram geometry.
For parallel beam this is just `ndgrid(sg.r, sg.ar)`
but for fan beam and mojette this involves more complicated computations.
"""
sino_geom_grid(sg::SinoPar) = ndgrid(sg.r, sg.ar)

function sino_geom_grid(sg::SinoFan)
    gamma = sg.gamma
    rad = sg.dso * sin.(gamma) + sg.source_offset * cos.(gamma)
    rg = repeat(rad, 1, sg.na) # [nb na]
    return (rg, gamma .+ sg.ar') # [nb na] phi = gamma + beta
end

function sino_geom_grid(sg::SinoMoj)
    phi = sg.ar
    # trick: ray_spacing aka ds comes from dx which is sg.d for mojette
    wb = (sg.nb - 1)/2 + sg.offset
    dt = sg.d_ang # [na]
    pos = ((0:(sg.nb-1)) .- wb) * dt' # [nb na]
    return (pos, repeat(phi', sg.nb, 1))
end

=#

"""
    show(io::IO, ::MIME"text/plain", sg::SinoGeom)
"""
function Base.show(io::IO, ::MIME"text/plain", sg::SinoGeom)
    println(io, "$(typeof(sg)) :")
    for f in propertynames(sg)
        p = getproperty(sg, f)
        t = typeof(p)
        println(io, " ", f, "::", t, " ", p)
    end
end

#sino_geom_how(sg::SinoGeom{G}) where {G} = G

#=
# Extended properties

sino_geom_fun0 = Dict([
    (:help, sg -> sino_geom_help()),

    (:dim, sg -> (sg.nb, sg.na)),
    (:w, sg -> (sg.nb-1)/2 + sg.offset),
    (:ones, sg -> ones(Float32, sg.dim)),
    (:zeros, sg -> zeros(Float32, sg.dim)),

    (:dr, sg -> ((sg isa SinoMoj) ? NaN : sg.d)),
    (:ds, sg -> sg.dr),
    (:r, sg -> sg.d * ((0:sg.nb-1) .- sg.w)),
    (:s, sg -> sg.r), # sample locations ('radial')

    (:gamma, sg -> sino_geom_gamma(sg)),
    (:gamma_max, sg -> maximum(abs.(sg.gamma))),
    (:orbit_short, sg -> 180 + 2 * rad2deg(sg.gamma_max)),
    (:ad, sg -> (0:sg.na-1)/sg.na * sg.orbit .+ sg.orbit_start),
    (:ar, sg -> deg2rad.(sg.ad)),

    (:rfov, sg -> sino_geom_rfov(sg)),
    (:xds, sg -> sino_geom_xds(sg)),
    (:yds, sg -> sino_geom_yds(sg)),
    (:dso, sg -> sg.dsd - sg.dod),
    (:grid, sg -> sino_geom_grid(sg)),
    (:plot_grid, sg -> ((plot::Function) -> sino_geom_plot_grid(sg, plot))),

    (:plot!, sg ->
        ((plot!::Function ; ig=nothing) -> sino_geom_plot!(sg, plot! ; ig=ig))),
    (:shape, sg -> (((x::AbstractArray) -> reshaper(x, sg.dim)))),
    (:taufun, sg -> ((x,y) -> sino_geom_taufun(sg,x,y))),
    (:unitv, sg -> ((;kwarg...) -> sino_geom_unitv(sg; kwarg...))),

    # angular dependent d for :moj
    (:d_moj, sg -> (ar -> sg.d * max(abs(cos(ar)), abs(sin(ar))))),
    (:d_ang, sg -> sg.d_moj.(sg.ar)),

    # functions that return new geometry:

    (:down, sg -> (down::Int -> downsample(sg, down))),
    (:over, sg -> (over::Int -> sino_geom_over(sg, over))),

    ])


# Tricky overloading here!

Base.getproperty(sg::SinoGeom, s::Symbol) =
        haskey(sino_geom_fun0, s) ? sino_geom_fun0[s](sg) :
        getfield(sg, s)

Base.propertynames(sg::SinoGeom) =
    (fieldnames(typeof(sg))..., keys(sino_geom_fun0)...)


"""
    reshaper(x::AbstractArray, dim:Dims)
Reshape `x` to size `dim` with `:` only if needed
"""
reshaper(x::AbstractArray, dim::Dims) =
    (length(x) == prod(dim)) ? reshape(x, dim) : reshape(x, dim..., :)


"""
    sino_geom_ge1()
sinogram geometry for GE lightspeed system
These numbers are published in IEEE T-MI Oct. 2006, p.1272-1283 wang:06:pwl
"""
function sino_geom_ge1( ;
    na::Int = 984,
    nb::Int = 888,
    orbit::Union{Symbol,Real} = 360,
    units::Symbol = :mm, # default units is mm
    kwarg...,
)

    if orbit === :short
        na = 642 # trick: reduce na for short scans
        orbit = na/984*360
    end

    scale = units === :mm ? 1 :
            units === :cm ? 10 :
            throw("units $units")

    return sino_geom_fan( ; units=units,
        nb=nb, na=na,
        d = 1.0239/scale, offset = 1.25,
        dsd = 949.075/scale,
        dod = 408.075/scale,
        dfs = 0, kwarg...,
    )
end
=#


"""
todo
    sino_geom_ge1()
sinogram geometry for GE lightspeed system
These numbers are published in IEEE T-MI Oct. 2006, p.1272-1283 wang:06:pwl
"""
function SinoFan(::Val{:ge1} ;
    na::Int = 984,
    nb::Int = 888,
    offset::Real = 1.25,
    orbit::Union{Symbol,Real} = 360,
    units::Symbol = :mm,
    dfs::RealU = 0,
    kwarg...,
)

    if orbit === :short
        na = 642 # trick: reduce na for short scans
        orbit = na/984*360
    end

    scale = units === :mm ? 1 :
            units === :cm ? 10 :
            throw("units $units")

    return SinoFanArc( ; nb, na, offset,
        d = 1.0239/scale,
        dsd = 949.075/scale,
        dod = 408.075/scale,
        kwarg...,
    )
end


function rays(sg::SinoPar)
end

# tests, todo

using Test

using Unitful: mm, °
d = 2mm
orbit = 180.0°
#d = 2; orbit = 180
sg = SinoMoj(; d, orbit)
sg = SinoFanArc(; d, orbit)
sg = SinoFanFlat(; d, orbit)
sg = SinoFan(Val(:ge1))

function _test(geo)
    @inferred geo(; d, orbit)
    args = 128, 100, d, orbit, 0*orbit, 0, 2d
    @inferred geo{Int,Float64}(args...)
    sg = @inferred geo(args...)
    downsample(sg, 2)
end

for geo in (SinoPar, SinoMoj)
    _test(geo)
end
