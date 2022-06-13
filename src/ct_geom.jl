#=
ct_geom.jl
CT geometry for 3D image reconstruction
2022-05-11, Jason Hu, Jeff Fessler translated from ct_geom.m in MIRT
=#
using PyPlot
using MIRT

export CtGeom
export CtFanPar, CtFan, CtFanArc, CtFanFlat
export ct_geom_help

abstract type CtGeom end
abstract type CtFanParallel end #line 238 of ct_geom
abstract type CtFan <: CtGeom end

# case 'par' for matlab
mutable struct CtFanPar <: CtFan
    units::Symbol             # :nothing | :mm etc.
    ns::Int                   #
    nt::Int
    na::Int
    orbit::Float64
    orbit_start::Float64
    ds::Float64
    dt::Float64
    offset_s::Float64
    offset_t::Float64
    dsd::Float64
    #dso::Float64 value is infinite
    dod::Float64
    down::Float64
    pitch::Float64
    user_source_zs
    source_z0::Float64
end

#case st.dfs==0 arc for matlab
mutable struct CtFanArc <: CtFan
    units::Symbol
    ns::Int
    nt::Int
    na::Int
    orbit::Float64
    orbit_start::Float64
    ds::Float64
    dt::Float64
    offset_s::Float64
    offset_t::Float64
    dsd::Float64
    #dso::Float64
    dod::Float64
    down::Float64
    pitch::Float64
    user_source_zs
    source_z0::Float64
    #not included: down, pitch, dfs, nframe, frame, source_z0, user_source_zs
end

#case isinf(st.dfs) flat for matlab
mutable struct CtFanFlat <: CtFan
    units::Symbol
    ns::Int
    nt::Int
    na::Int
    orbit::Float64
    orbit_start::Float64
    ds::Float64
    dt::Float64
    offset_s::Float64
    offset_t::Float64
    dsd::Float64
    #dso::Float64
    dod::Float64
    down::Float64
    pitch::Float64
    user_source_zs
    source_z0::Float64
end

#help function
function ct_geom_help( ; io::IO = isinteractive() ? stdout : devnull )
    print(io,
    "\n
    Derived values

    st.s			s sample locations
	st.t			t sample locations
	st.gamma		[nb] gamma sample values [radians]
	st.gamma_s		gamma values given s values
	st.gamma_max		half of fan angle [radians] (if offset_s=0)
	st.gamma_max_abs	half of fan angle [radians] (recommended)
	st.ws			(ns-1)/2 + st.offset_s
	st.wt			(nt-1)/2 + st.offset_t
	st.ad			[na] source angles in degrees
	st.ar			[na] source angles in radians
	st.dim			dimensions: [st.ns st.nt st.na]
	st.downsample(down)	reduce sampling by integer factor
	st.ones			ones(ns,nt,na, 'single')
	st.zeros		zeros(ns,nt,na, 'single')
	st.rad			[ns] radial distance of each ray from origin
	st.rmax			max radius within FOV
	st.footprint_size(ig)	max footprint width in 's'
	st.cone_angle		(half) cone angle on axis (s=0): +/- angle
	st.zfov			axial FOV
	st.source_zs		[na] z-locations of source for each view

    Methods

    st.down(down)   reduce sampling by integer factor
    st.plot2(ig)    plot the 2D geometry of the image
    st.plot3(ig)    plot the 3D geometry of the cone beam setup
    st.unitv(is,it,ia)  unit vector with single nonzero element

    ")
end

"""
Constructor for CtGeom
Create the CT geometry that describes the sampling characteristics of a
CT scan using the parallel, fan beam arc, or flat fan beam geometry
Use sino_geom() instead for the 2D version
"""
function ct_geom(how::Symbol ; kwarg...)
    #line 177 of sino_geom.jl? ge1 type
    if how === :par
        sg = ct_geom_par( ; kwarg...)
    elseif how === :flat
        sg = ct_geom_flat( ; kwarg...)
    elseif how === :arc
        sg = ct_geom_arc( ; kwarg...)
    elseif how === :ge1
        sg = ct_geom_ge1( ; kwarg...)
    else
        throw("unknown ct type $how")
    end
    return sg
end

#line 193 of sino_geom.jl?

#functions that are common to all
function ct_geom_orbit_short(st::CtGeom)
    return 180 + 2 * rad2deg(st.gamma_max)
end

function ct_geom_dim(st::CtGeom)
    return [st.ns, st.nt, st.na]
end

function ct_geom_flat( ;
    units::Symbol = :none,
    ns::Int = 128,
    nt::Int = 64,
    na::Int = 64,
    orbit::Union{Symbol,Real} = 360,
    orbit_start::Real = 0,
    ds::Real = 1,
    dt::Real = 1,
    offset_s::Real = 0,
    offset_t::Real = 0,
    dsd::Float64 = 4.0*ns,
    dod::Float64 = 0.0,
    down::Int = 1,
    pitch::Float64 = 0.0,
    user_source_zs = [],
    source_z0::Real = 0
    )
    st = CtFanFlat(units, ns, nt, na, orbit, orbit_start, ds, dt, offset_s, offset_t, dsd, dod, down, pitch, user_source_zs, source_z0)
    if down != 1
        return st.downsample(down)
    end
    return st # need to downsample this later
end

function ct_geom_arc( ;
    units::Symbol = :none,
    ns::Int = 128,
    nt::Int = 64,
    na::Int = 64,
    orbit::Union{Symbol,Real} = 360,
    orbit_start::Real = 0,
    ds::Real = 1,
    dt::Real = 1,
    offset_s::Real = 0,
    offset_t::Real = 0,
    dsd::Float64 = 4.0*ns,
    dod::Float64 = 0.0,
    down::Int = 1,
    pitch::Float64 = 0.0,
    user_source_zs = [],
    source_z0::Real = 0
    )
    st = CtFanArc(units, ns, nt, na, orbit, orbit_start, ds, dt, offset_s, offset_t, dsd, dod, down, pitch, user_source_zs, source_z0)
    return st.downsample(down)
end

function ct_geom_par( ;kwarg...)
    error("Not implemented yet")
    return 0
end

function ct_geom_ad(st::CtGeom)
    return [0:st.na-1;] / st.na * st.orbit .+ st.orbit_start
end

function ct_geom_ar(st::CtGeom)
    return deg2rad.(ct_geom_ad(st))
end

function ct_geom_shape(st::CtGeom, sino)
    error("Not implemented yet")
end

function ct_geom_rmax(st::CtFanPar)
    #line 262 of ct_geom.m doesn't make sense
    #check if st.type is fan, but then if it's parallel?
    smax = maximum(abs.(st.s))
    return smax
end

function ct_geom_rmax(st::CtFanArc)
    gamma_max_abs = ct_geom_gamma_max_abs(st)
    return st.dso * sin(gamma_max_abs)
end

function ct_geom_rmax(st::CtFanFlat)
    gamma_max_abs = ct_geom_gamma_max_abs(st)
    return st.dso * sin(gamma_max_abs)
end

function ct_geom_xds(st::CtFanArc)
    ss = st.s
    return st.dsd * sin.(st.gamma_s(ss))
end

function ct_geom_xds(st::CtFanFlat)
    return st.s
end

function ct_geom_xds(st::CtFanParallel)
    return st.s
end

function ct_geom_yds(st::CtFanArc)
    ss = st.s
    return st.dso .- st.dsd * cos.(st.gamma_s(ss))
end

function ct_geom_yds(st::CtFanFlat)
    ss = st.s
    return -st.dod * ones(size(ss), class(ss))
end

function ct_geom_yds(st::CtFanParallel)
    ss = st.s
    return zeros(size(ss), class(ss))
end

function ct_geom_cone_angle(st::CtFanParallel)
    return 0
end

function ct_geom_cone_angle(st::CtFanArc)
    return atan((st.nt * st.dt)/2 / st.dsd)
end

function ct_geom_cone_angle(st::CtFanFlat)
    return atan((st.nt * st.dt)/2 / st.dsd)
end

function ct_geom_zfov(st::CtFanParallel)
    return st.nt * st.dt
end

function ct_geom_cone_angle(st::CtFanArc)
    return st.dso / st.dsd * st.nt * st.dt
end

function ct_geom_cone_angle(st::CtFanFlat)
    return st.dso / st.dsd * st.nt * st.dt
end

function ct_geom_source_dz_per_view(st::CtGeom)
    if st.na == 1
        return 0
    end
    if st.pitch == 0
        return 0
    end
    if length(st.orbit) != 1 || st.orbit == 0
        error("ERROR OCCURRED")
        return 0
    end
    na_per_360 = st.na * (360 / st.orbit)
    out = st.pitch * st.zfov / na_per_360
end

function ct_geom_source_zs(st::CtGeom)
    if !isempty(st.user_source_zs)
        return st.user_source_zs
    else
        source_dz = st.source_dz_per_view
        return [0:st.na-1;] * source_dz .+ st.source_z0
    end
end

function ct_geom_gamma_s(st::CtFanArc, ss)
    return ss / st.dsd
end

function ct_geom_gamma_s(st::CtFanFlat, ss)
    return atan.(ss / st.dsd)
end

function ct_geom_gamma(st::CtGeom)
    return ct_geom_gamma_s(st, st.s)
end

function ct_geom_gamma_max(st::CtGeom)
    return maximum(st.gamma)
end

function ct_geom_gamma_max_abs(st::CtGeom)
    return maximum(abs.(st.gamma))
end

function ct_geom_footprint_size(st::CtFanPar, ig)
    di = sqrt(ig.dx^2 + ig.dy^2)
    smax = maximum(abs.(st.s))
    rfov = max(ig.nx * ig.dx, ig.ny * ig.dy) / 2
    dso = st.dso
    dsd = st.dsd
    if rfov > 0.99 * dso
        error("BAD DSO ERROR")
        return 1 / 0
    end
    return di / st.ds
end

function ct_geom_footprint_size(st::CtFanFlat, ig)
    di = sqrt(ig.dx^2 + ig.dy^2)
    smax = maximum(abs.(st.s))
    rfov = max(ig.nx * ig.dx, ig.ny * ig.dy) / 2
    dso = st.dso
    dsd = st.dsd
    if rfov > 0.99 * dso
        error("BAD DSO ERROR")
        return 1 / 0
    end
    return di / st.ds * sqrt(dsd^2 + smax^2) / (dso-rfov)
end

function ct_geom_footprint_size(st::CtFanArc, ig)
    di = sqrt(ig.dx^2 + ig.dy^2)
    smax = maximum(abs.(st.s))
    rfov = max(ig.nx * ig.dx, ig.ny * ig.dy) / 2
    dso = st.dso
    dsd = st.dsd
    if rfov > 0.99 * dso
        error("BAD DSO ERROR")
        return 1 / 0
    end
    return di / st.ds * dsd / (dso-rfov)
end

function ct_geom_downsample(si, down)
    st = si
    #st.down = st.down .* down
    down = st.down #I only have this to give the right results
    if length(down) == 1
        down_s = down
        down_t = down
        down_a = down
    elseif length(down) == 3
        down_s = down[1]
        down_t = down[2]
        down_a = down[3]
    else
        print("Error in down elements")
        return 1 / 0
    end

    st.ns = 4 * ceil(st.ns / down_s / 4)
    st.nt = 2 * ceil(st.nt / down_t / 2)

    user_spec = false
    if !isempty(st.user_source_zs)
        st.user_source_zs = st.user_source_zs[1:down_a:st.na]
        user_spec = true
    end

    if length(st.orbit_start) > 1
        st.orbit_start = st.orbit_start[1:down_a:si.na]
        user_spec = true
    end

    st.ds = st.ds * down_s
    st.dt = st.dt * down_t

    if user_spec
        st.na = length([1:down_a:st.na])
    elseif all(diff(si.source_zs) == 0)
        st.na = max(1, round(si.na / down_a))
    else
        nturn = si.orbit / 360
        na1 = si.na / nturn

        tol = 0.0001
        if abs(round(na1)-na1) < tol
            na2 = round(na1 / down_a)
            tmp = nturn * na2
            if abs(round(tmp)-tmp) < tol
                st.na = round(tmp)
            else
                st.na = floor(tmp)
                st.orbit = 360 / na2 * st.na
            end
        else
            st.na = round(si.na / down_a)
        end
    end
    return st
end

function ct_geom_unitv(st, is, it, ia)
    out = st.zeros
    if isempty(is)
        is = floor(st.ns/2 + 1)
        it = floor(st.nt/2 + 1)
        ia = 1
    end
    out[is, it, ia] = 1
    return out
end

#just to make the other parts of the code more convenient that use dfs
function ct_geom_getDFS(st::CtGeom)
    if isa(st, CtFanArc)
        return 0
    elseif isa(st, CtFanFlat)
        return Inf
    elseif isa(st, ct_geom2)
        return st.dfs
    else
        error("DFS Not implemented yet")
    end
    return 0
end

#plot 2D geometry of the setup
function ct_geom_plot2(st::CtGeom, ig::ImageGeom)
    t = LinRange(0, 2*pi, 1001)

    if isa(st, CtFanParallel)
        error("Not implemented yet")
    end
    beta0 = deg2rad(st.orbit_start)
    rot = [cos(beta0) sin(beta0); -sin(beta0) cos(beta0)]
    p0 = st.dso * [-sin(beta0); cos(beta0)]
    pd = rot' * [st.xds'; st.yds']
    rfov = st.dso * sin(maximum(abs.(st.gamma)))

    Plots.scatter([p0[1]],[p0[2]], linestyle = :solid)
    Plots.plot!(st.dso * cos.(t), st.dso * sin.(t), linecolor = :red)
    Plots.plot!(pd[1,:], pd[2,:], linestyle = :solid, linecolor = :yellow)

    xmin = minimum(ig.x)
    xmax = maximum(ig.x)
    ymin = minimum(ig.y)
    ymax = maximum(ig.y)
    Plots.plot!([xmax, xmin, xmin, xmax, xmax], [ymax, ymax, ymin, ymin, ymax], linecolor = :green)

    Plots.plot!([pd[1,1], p0[1], last(pd[1,:])], [pd[2,1], p0[2], last(pd[2,:])], linestyle = :dash)
    Plots.plot!(rfov*cos.(t), rfov*sin.(t), linecolor = :magenta, linestyle = :dot)
end

#helper functions for ct_geom_plot3
function r100(x)
    return 100 * ceil(maximum(abs.(x))/100)
end

function r10(x)
    return 10 * ceil(maximum(abs.(x))/10)
end

#3D helical geometry plot, saves plot to myfig.png
function ct_geom_plot3(st::CtGeom, ig::ImageGeom)
    if isinf(st.dso) || isinf(st.dsd)
        error("Parallel beam Not done")
    end

    t1 = -st.dso * sin.(st.ar)
    t2 = st.dso * cos.(st.ar)
    t3 = st.source_zs

    fig = figure()
    ax = Axes3D(fig)
    ax.plot3D(t1, t2, t3, color = "yellow", linestyle = "dotted")

    ax.plot3D([-1, 1] * r100(t1), [0,0], [0,0])
    ax.plot3D([0,0], [-1, 1]*r100(t2), [0,0])
    ax.plot3D([0,0], [0,0], [-1,1]*r10(t3))

    #skip the part about ztick line 740
    xmin = minimum(ig.x)
    xmax = maximum(ig.x)
    ymin = minimum(ig.y)
    ymax = maximum(ig.y)
    zmin = minimum(ig.z)
    zmax = maximum(ig.z)
    ax.plot3D([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], [zmin, zmin, zmin, zmin, zmin], color = "green", linestyle = "dashed")
    ax.plot3D([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], [zmax, zmax, zmax, zmax, zmax], color = "green", linestyle = "dashed")
    ax.plot3D(xmin*[1,1], ymin*[1,1], [zmin, zmax], color = "green", linestyle = "dashed")
    ax.plot3D(xmax*[1,1], ymin*[1,1], [zmin, zmax], color = "green", linestyle = "dashed")
    ax.plot3D(xmin*[1,1], ymax*[1,1], [zmin, zmax], color = "green", linestyle = "dashed")
    ax.plot3D(xmax*[1,1], ymax*[1,1], [zmin, zmax], color = "green", linestyle = "dashed")
    ax.plot3D(xmin*ones(ig.nz), ymin*ones(ig.nz), ig.z, color = "green", linestyle = "dotted")

    detcolor = ["red", "cyan", "magenta"]
    if st.na == 1
        ia_list = 1
    elseif st.na == 2
        ia_list = [1, st.na]
    else
        ia_list = [1, 1+floor(Int, st.na/2), st.na]
    end

    for ia = ia_list
        src = [t1[ia], t2[ia], t3[ia]]
        unit_vec = [-src[1], -src[2], 0]/sqrt(src[1]^2 + src[2]^2)
        det_cen_loc = src + st.dsd * unit_vec
        #ax.plot3D(det_cen_loc[1], det_cen_loc[2], det_cen_loc[3], )

        rot = st.ar[ia]
        rot = [cos(rot) -sin(rot); sin(rot) cos(rot)]
        pd = rot * [st.xds'; st.yds']

        for it = 1:st.nt
            indices = findall(x -> x == ia, ia_list)
            indices = indices[1]
            ax.plot3D(pd[1,:], pd[2,:], src[3] .+ st.t[it] * ones(st.ns), color = detcolor[indices])
        end
    end
    ax.view_init(22, -200)
    PyPlot.savefig("myfig.png")
end

ct_geom_fun0 = Dict([
    (:help, st -> ct_geom_help()),

    (:dim, st -> (st.ns, st.nt, st.na)),
    (:ws, st -> (st.ns-1)/2 + st.offset_s),
    (:wt, st -> (st.nt-1)/2 + st.offset_t),
    (:ones, st -> ones(Float64, st.dim)),
    (:zeros, st -> zeros(Float64, st.dim)),

    (:s, st -> st.ds * ([0:st.ns-1;] .- st.ws)),
    (:t, st -> st.dt * ([0:st.nt-1;] .- st.wt)),
    (:dso, st -> st.dsd - st.dod),
    (:dfs, st -> ct_geom_getDFS(st)),

    (:gamma, st -> ct_geom_gamma(st)),
    (:gamma_s, st -> (ss -> ct_geom_gamma_s(st,ss))),
    (:gamma_max, st -> ct_geom_gamma_max(st)),
    (:gamma_max_abs, st -> ct_geom_gamma_max_abs(st)),
    (:xds, st -> ct_geom_xds(st)),
    (:yds, st -> ct_geom_yds(st)),

    (:cone_angle, st -> ct_geom_cone_angle(st)),
    (:zfov, st -> ct_geom_zfov(st)),
    (:source_dz_per_view, st -> ct_geom_source_dz_per_view(st)),
    (:source_zs, st -> ct_geom_source_zs(st)),
    (:orbit_short, st -> ct_geom_orbit_short(st)),
    (:downsample, st -> (down::Int -> ct_geom_downsample(st, down))),

    (:ad, st -> ct_geom_ad(st)),
    (:ar, st -> ct_geom_ar(st)),
    (:rad, st -> ct_geom_rad(st)),
    (:rmax, st -> ct_geom_rmax(st)),
    (:shape, st -> ct_geom_shape(st)),
    (:footprint_size, st -> ct_geom_footprint_size(st)),

    (:plot_grid, st -> ((plot::Function) -> ct_geom_plot_grid(st, plot))),
    (:plot2, st -> ((ig::ImageGeom) -> ct_geom_plot2(st, ig))),
    (:plot3, st -> ((ig::ImageGeom) -> ct_geom_plot3(st, ig)))
])

Base.getproperty(st::CtGeom, s::Symbol) =
    haskey(ct_geom_fun0, s) ? ct_geom_fun0[s](st) :
    getfield(st, s)

Base.propertynames(st::CtGeom) =
    (fieldnames(typeof(st))..., keys(ct_geom_fun0)...)

# CT geometry for GE lightspeed system
# These numbers are published in IEEE T-MI Oct. 2006, p.1272-1283 wang:06:pwl
function ct_geom_ge1( ;
    na::Int = 1)
    st = ct_geom(:arc ; ns = 888, nt = 64, na = 984, ds = 1.0239, dt = 1.0964, offset_s = 1.25, dsd = 949.075, dod = 408.075)
    return st
end
