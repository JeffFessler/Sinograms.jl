# fbp2.jl

export fbp2

using ImageGeoms: ImageGeom
using Sinograms: SinoGeom, SinoPar, SinoFan, SinoMoj

abstract type FBPplan end

struct Moj{Th,Tg}
    H::Th
    G::Tg
end

Moj() = Moj{AbstractMatrix{<:Real},AbstractMatrix{<:Real}}(zeros(Real,1,1),zeros(Real,1,1))
#Moj(a,b) = Moj{AbstractMatrix{<:Real},AbstractMatrix{<:Real}}(zeros(Real,1,1),zeros(Real,1,1)) # todo method overwritten warning


#struct NormalPlan{S, I, W, P} <: FBPplan
struct NormalPlan{
    S <: SinoGeom,
    I <: ImageGeom,
    W <: Window,
    P <: AbstractMatrix{<:Real},
} <: FBPplan
    sg::S
    ig::I
    window::W
    parker_weight::P
#   moj::Moj # todo

#=
    function NormalPlan(
        sg::S,
        ig::I,
        window::W,
        parker_weight::P,
    ) where {S <: SinoGeom, I <: ImageGeom,
        W <: Window, P <: AbstractMatrix{<:Real}}
        return NormalPlan{S,I,W,P}(sg, ig, window, parker_weight) #, Moj())
    end
=#
end


struct DfPlan <: FBPplan
    sg::SinoGeom
    ig::ImageGeom
    window::Window
    parallel_beam_parker_weight::AbstractMatrix{<:Real}
end

struct MojPlan <: FBPplan
    sg::SinoGeom
    ig::ImageGeom
    window::Window
    parallel_beam_parker_weight::AbstractMatrix{<:Real}
    moj::Moj
end

struct TabPlan <: FBPplan
    sg::SinoGeom
    ig::ImageGeom
    window::Window
    parallel_beam_parker_weight::AbstractMatrix{<:Real}
end

reale = (x) -> (@assert x ≈ real(x); real(x))



"""
    plan = fbp2(sg, ig; how=:normal, window=Window())

FBP 2D tomographic image reconstruction for parallel-beam or fan-beam cases,
with either flat or arc detector for fan-beam case.

To use this, you first call it with the sinogram and image geometry.
The routine returns the initialized "plan". Thereafter, to
to perform FBP reconstruction, you call this routine with the plan
(perhaps numerous times for the same geometry).


in
- `sg::SinoGeom`
- `ig::ImageGeom`

options
- `how::Symbol`             how to reconstruct
    * `:normal`             default
    * `:mojette`            use mojette rebinning and Gtomo2_table
- `window::Window` e.g. `Window(Hamming(0.5))`; default `Window()`
-`T::DataType`              type of sino elements (default: `Float32`)

out
- `plan::FBPplan`            initialized plan

"""
function fbp2(
    sg::SinoGeom,
    ig::ImageGeom;
    how::Symbol = :normal,
    window::Window = Window(),
    T::DataType = Float32,
#   nthread::Int = Threads.nthreads()
)

    if how === :normal
        plan = fbp_setup_normal(sg, ig, window, T)
#   elseif how === :dsc
#       plan = fbp_setup_dsc(sg, ig, how, window)
#   elseif how === :df
#       plan = fbp_setup_df(sg, ig, how, window)
#   elseif how === :mojette
#       plan = fbp_setup_moj(sg, ig, how, window)
#   elseif how === :table
#       plan = fbp_setup_tab(sg, ig, how, window)
    else
        throw("unknown type: $how")
    end

    return plan
end

function fbp2_par_parker_wt(sg::SinoGeom)
    orbit = abs(sg.orbit)
    na = sg.na
    nb = sg.nb
    ad = abs(sg.ad - sg.orbit_start)

    if orbit < 180
        @warn("orbit $orbit < 180")
        return sg.ones
    end

    orbit > 360 && throw("only 180 ≤ orbit ≤ 360 supported for Parker weighting")
    extra = orbit - 180 # extra beyond 180

    wt = ones(T,na)

    ii = ad .< extra
    wt[ii] = abs2.(sin(ad[ii] ./ extra .* pi ./ 2))
    ii = ad .>= 180
    wt[ii] = abs2.(sin((orbit .- ad[ii]) ./ extra .* pi ./ 2))
    wt *= orbit / 180 # trick because of the back-projector normalization
    return repeat(wt, nb, 1) # [nb na] sinogram sized
end


function fbp_setup_normal(
    sg::SinoPar,
    ig::ImageGeom,
    window::Window,
    T::DataType,
)
    weight = ones(T, sg.nb, sg.na)
    if abs(sg.orbit) != 180 && abs(sg.orbit) != 360
        weight = fbp2_par_parker_wt(sg)
    end
    return NormalPlan(sg, ig, window, weight)
end


function fbp_setup_normal(
    sg::SinoFan,
    ig::ImageGeom,
    window::Window,
    T::DataType,
)
    weight = ones(T, sg.nb, sg.na)
    sg.orbit != 360 && @warn("short-scan fan-beam Parker weighting not done")
    return NormalPlan(sg, ig, window, weight)
end


function fbp_setup_normal(
    sg::SinoMoj,
    ig::ImageGeom,
    window::Window,
    T::DataType,
)
    weight = ones(T, sg.nb, sg.na)
    if sg.dx == abs(ig.dx)
        throw("not done")
#       plan.moj.G = Gtomo2_table(sg, ig, ["mojette,back1"], nthread=nthread)
    else
        d = sg.dx
        dx = ig.dx
        @warn("mojette sinogram with d=$d vs image with dx=$dx")
    end

    plan.moj.H = fbp2_make_sino_filter_moj(sg.nb, sg.na, sg.dx, sg.orbit, sg.orbit_start, window)
    return NormalPlan{SinoMoj}(sg, ig, window, weight, moj)
end


"""
    image, sino_filt = fbp2(plan, sino)

Filtered back-projection (FBP) reconstruction,
returning image and filtered sinogram

in
- `plan::FBPplan`
- `sino::AbstractMatrix{<:Number}`

out
- `image::Matrix{<:Number}`       reconstructed image(s)
- `sino_filt::Matrix{<:Number}`   filtered sinogram(s)

"""
function fbp2(plan::NormalPlan, sino::AbstractMatrix{<:Number})

    plan.sg.dim != size(sino) && throw("bad sino size")

    if plan.sg isa SinoPar
        sino .*= plan.parker_weight

        sino,_,_,_ = fbp_sino_filter(plan.sg, sino, window = plan.window)

        #image = fbp_back(plan.sg, plan.ig, sino, do_r_mask=true)
        image = fbp_back(plan.sg, plan.ig, sino)

    elseif plan.sg isa SinoFan

        dfs = plan.sg.dfs

        dfs != 0 && ~isinf(dfs) && throw("only arc or flat fan done")

        sino .*= fbp_sino_weight(plan.sg)

        sino,_,_,_ = fbp_sino_filter(
            plan.sg, sino,
            window = plan.window,
        )
        image = fbp_back(plan.sg, plan.ig, sino)

    elseif plan.sg isa SinoMoj #TODO (incomplete)

        sino = fbp2_apply_sino_filter_moj(sino, plan.moj.H)

        if plan.sg.dx == abs(plan.ig.dx)
            image = plan.moj.G' * sino # backproject
            image = image * (pi / plan.sg.na) # account for "dphi" in integral
        else # revert to conventional pixel driven
            ig = plan.ig
            sg = plan.sg
            arg1 = [uint8(ig.mask), ig.dx, ig.dy, ig.offset_x, sign(ig.dy) * ig.offset_y] # trick: old backproject
            arg2 = [sg.d(1:sg.na), sg.offset, sg.orbit, sg.orbit_start]
# todo      image = jf_mex("back2", arg1[:], arg2[:], int32(arg.nthread), single(sino))
            image = image .* plan.ig.mask
        end

    else
        throw("not done")
    end

    return image, sino
end


# fbp2_recon_moj
function fpb2(plan::MojPlan, sino::AbstractMatrix{<:Number})
    moj = plan.moj
    # bin sinogram to mojette sampling, whether it is fan or parallel
    msino = moj.ob_rebin * sino
    msino = fbp2_apply_sino_filter_moj(msino, moj.H) # filter mojette sinogram
    image = moj.G' * msino # backproject
    return image * (π / size(sino,2)) # account for "dphi" in integral
end

# todo: needs broadcasts?
function fbp_make_sino_filter_moj(nb, na, dx, orbit, orbit_start, window)
    ang = deg2rad(orbit_start .+ (0:na-1)./na .* orbit)
    npad = 2^ceil(log2(2*nb-1)) # padded size

    dr = dx * max(abs(cos(ang)), abs(sin(ang)))
    if true
        junk, H = fbp2_sino_filter(:flat, ones(nb,1), ds=1, window=:window, decon1=0)

        # trick: ramp filter usually has 1/dr^2 in it, but convolution by sum
        # requires "dr" so we need dr / dr^2 = 1 / dr
        return H * (1 ./ dr) # [npad,na]

    else # trick: make sure cutoff frequencies match
        u0 = 1/2/dx # standard cutoff for coarsest sampling
        r = -(npad/2):(npad/2-1) * dr # [nb na]
        h = u0^2 .* (2 .* sinc.(2*u0*r) - sinc.(u0*r).^2)
        h = h .* repeat(dr, [npad 1]) # extra dr for discrete-space convolution
   #    clf, plot(r, h, '.'), keyboard
        return reale(fft(fftshift(h,1), 1))
   #    H = fbp_apodize(H, nb, window)
        !isempty(window) && throw("window not done yet due to dr")
    end
end


function fbp2_apply_sino_filter_moj(sino, H)
    nb = size(sino,1)
    npad = 2^ceil(log2(2*nb-1)) # padded size
    sinopad = [sino; zeros(npad-nb,size(sino,2))] # padded sinogram
    sino = ifft(reale(fft(sinopad, 1) .* H), 1)
    sino = sino[1:nb,:]
end
