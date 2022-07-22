# fbp2/ramp.jl

export fbp_ramp, ramp_flat, ramp_arc

#using Sinograms: SinoFanFlat, SinoFanArc, SinoPar, RealU


"""
    h, n = fbp_ramp(sg::SinoGeom, N::Int)
'ramp-like' filters for parallel-beam and fan-beam FBP reconstruction.
This sampled band-limited approach avoids the aliasing that would be
caused by sampling the ramp directly in the frequency domain.

in
- `sg::SinoGeom`
- `N::Int` : # of samples (must be even)

out
- `h::Vector{<:Number}` : samples of band-limited ramp filter
- `n::UnitRange{Int64}` : -(N÷2):(N÷2-1)

"""
fbp_ramp(sg::SinoPar, N::Int) = ramp_flat(N, sg.d)
fbp_ramp(sg::SinoFanFlat, N::Int) = ramp_flat(N, sg.d)
fbp_ramp(sg::SinoFanArc, N::Int) = ramp_arc(N, sg.d, sg.dsd)


function _ramp_arc(n::Int, ds::RealU, dsd::RealU)
    R = promote_type(eltype(ds), eltype(dsd), eltype(1f0 * oneunit(ds))) # at least Float32
    R = eltype(1 / oneunit(R)^2)
    h = n == 0 ? 0.25 / abs2(ds) :
        isodd(n) ? -1 / abs2(π * dsd * sin(n * ds / dsd)) :
        zero(1 / abs2(ds))
    return R(h)
end


function _ramp_flat(n::Int, ds::RealU)
    R = promote_type(eltype(ds), eltype(1f0 * oneunit(ds))) # at least Float32
    R = eltype(1 / oneunit(R)^2)
    h = n == 0 ? 0.25 / abs2(ds) :
        isodd(n) ? -1 / abs2(π * n * ds) :
        zero(1 / abs2(ds))
    return R(h)
end


function ramp_arc(N::Int, ds::RealU, dsd::RealU)
    isodd(N) && throw("N must be even")

    (N/2 * ds / dsd > 0.9 * π/2) &&
        throw("angle is $(rad2deg(N/2 * ds / dsd)) degrees: too large;
        physically impossible arc geometry")

    n = -(N÷2):(N÷2-1)
    h = _ramp_arc.(n, ds, dsd)
    return h, n
end


function ramp_flat(N::Int, ds::RealU)
    isodd(N) && throw("N must be even")
    n = -(N÷2):(N÷2-1)
    h = _ramp_flat.(n, ds)
    return h, n
end
