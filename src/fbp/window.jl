# fbp/window.jl
# Windows for apodizing the ramp filter.

export fbp_window, Window
export NoWindow, Boxcar, Hamming, Hann, WindowVect

using FFTW: fftshift


abstract type AbstractWindowShape end

"""
    Window{S,T}
Data type for FBP apodizing windows.
"""
struct Window{S,T}
    shape::S
    cutoff::T
    function Window(
        shape::S = NoWindow(),
        cutoff::T = 1,
    ) where {S <: AbstractWindowShape, T <: Real}
        return new{S,T}(shape, cutoff)
    end
end



"""
    WindowVect{V} <: AbstractWindowShape
A user-specified window vector, constructed via `WindowVect(v::V)`,
where `v` is a `AbstractVector`.
Caution: `length(v)` must be appropriate for the padded sinogram size.
"""
struct WindowVect{V} <: AbstractWindowShape
    v::V
    function WindowVect(v::V) where {V <: AbstractVector{<:Real}}
         new{V}(v)
    end
end


"""
    fbp_window(w::Window, N::Int ; T = Float32)

Create an apodizing window of length `N` and `fftshift` it.
"""
function fbp_window(w::Window, N::Int ; T = Float32)
    width = w.cutoff * N
    win = [T(window(w, width, n)) for n in (0:N-1) .- N÷2]
    return fftshift(win)::Vector{T}
end

function fbp_window(
    w::Window{S},
    N::Int ;
    T::DataType = eltype(w.shape.v),
) where {S <: WindowVect}
    N == length(w.shape.v) || throw("N mismatch")
    return fftshift(T.(w.shape.v))
end


# specific window shapes

struct NoWindow <: AbstractWindowShape end
window(w::Window{NoWindow}, width::Real, n::Int) = 1

struct Boxcar <: AbstractWindowShape end
window(w::Window{Boxcar}, width::Real, n::Int) = abs(n) < (width / 2)

struct Hamming <: AbstractWindowShape end
window(w::Window{Hamming}, width::Real, n::Int) =
    (0.54 + 0.46 * cos(2π*n/width)) * (abs(n) < (width / 2))

struct Hann <: AbstractWindowShape end
window(w::Window{Hann}, width::Real, n::Int) =
    0.5 * (1 + cos(2π*n/width)) * (abs(n) < (width / 2))
