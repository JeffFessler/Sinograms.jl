#---------------------------------------------------------
# # [Tomography](@id 01-tomography)
#---------------------------------------------------------

#=
This page gives an overview of the Julia package
[`Sinograms.jl`](https://github.com/JeffFessler/Sinograms.jl).

This page was generated from a single Julia file:
[01-tomography.jl](@__REPO_ROOT_URL__/01-tomography.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`01-tomography.ipynb`](@__NBVIEWER_ROOT_URL__/01-tomography.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`01-tomography.ipynb`](@__BINDER_ROOT_URL__/01-tomography.ipynb).


# ### Setup

# Packages needed here.

using Sinograms: fbp
using MIRTjim: jim, prompt
using InteractiveUtils: versioninfo


# The following line is helpful when running this example.jl file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


#=
### Tomography

[Tomography](https://en.wikipedia.org/wiki/Tomography)
is the process of imaging cross sections of an object
without actually slicing the object.
There are many forms of tomography;
the description here
focuses on
[X-ray computed tomography (CT scans)](https://en.wikipedia.org/wiki/CT_scan).
(See
[SPECTrecon.jl](https://github.com/JeffFessler/SPECTrecon.jl)
for a Julia package related to
[SPECT](https://en.wikipedia.org/wiki/Single-photon_emission_computed_tomography)
imaging.)

In an
[X-ray CT imaging system](http://en.wikipedia.org/wiki/File:Ct-internals.jpg),
X-rays emitted from an X-ray source
are transmitted through an object
(e.g., a patient in medical CT)
towards a detector array.
The source and detector
[rotate rapidly](https://www.youtube.com/watch?v=2CWpZKuy-NE)
around the object.
The signal intensities recorded by the detector
are related
to the sizes and densities
of the object materials
between the source and each detector element.


## Radon transform

The mathematical foundation
for 2D X-ray CT imaging
is the
[Radon transform](https://en.wikipedia.org/wiki/Radon_transform),
an
[integral transform](https://en.wikipedia.org/wiki/Integral_transform)
that maps a 2D function
``f(x,y)``
into the collection of line integrals
through that function.
Here we describe each line
by its angle ``\phi``
measured counter-clockwise from the ``y`` axis,
and by the radial distance ``r`` of the line from the origin.
The collection of integrals
is called a sinogram.

Mathematically,
the Radon transform ``p(r,\phi)`` of ``f(x,y)`` is defined by
```math
p(r,\phi) = \int_{-\infty}^{\infty}
f(r \cos \phi - l \sin \phi,
r \sin \phi + l \cos \phi) \, \mathrm{d} l.
```

The Radon transform of the object that is 1 inside the unit disk
and 0 elsewhere is given by
```math
p_1(r,\phi) = 2 \sqrt{1 - r^2}, \ \mathrm{ for } \ |r| < 1.
```

By the Radon transform's translation and scaling properties,
the Radon transform
of a disk of radius ``r_0``
centered at ``(x_0,y_0)``
is given by
```math
p(r,\phi) = r_0 \, p_1
\left( \frac{r - (x_0 \cos \phi + y_0 \sin \phi)}{r_0}, ϕ \right).
```

## Sinogram example

Here is a display of that Radon transform
for a disk of radius ``3``
centered at coordinate ``(5,1)``.
Note that maximum value is approximately ``6``,
the length of the longest chord
through a disk of radius ``3``.
=#

nr = 128
dr = 20 / nr
r = ((-(nr-1)/2):((nr-1)/2)) * dr # radial sample locations
na = 130
ϕ = (0:(na-1))/na * π # angular samples
proj1 = r -> abs(r) < 1 ? 2 * sqrt(1 - r^2) : 0.
proj2 = (r, ϕ, x, y, r0) -> r0 * proj1(r/r0 - (x * cos(ϕ) + y * sin(ϕ))/r0)
sino = proj2.(r, ϕ', 5, 1, 3)
jimsino = (sino, title) -> jim(
    r, ϕ, sino; title, aspect_ratio=:none,
    xlabel = "r", ylabel = "ϕ", ylims=(0,π), yticks=([0, π], ["0", "π"]),
    yflip=false, xticks = [-10, 0, 2, 5, 8, 10],
    clim = (0, 6),
)
jimsino(sino, "Sinogram for one disk")


#=
As this figure illustrates,
the Radon transform of a unit disk
has a somewhat sinusoidal shape.
Indeed every point in the ``(x,y)`` plane
traces out a distinct sinusoid
in the sinogram.

The mapping from the object ``f(x,y)``
to data like the sinogram ``p(r,\phi)``
is called the "forward problem."


## Image reconstruction

[Image reconstruction](http://en.wikipedia.org/wiki/Image_reconstruction)
is the process of solving the
[inverse problem](https://en.wikipedia.org/wiki/Inverse_problem)
of recovering an estimate ``\hat{f}(x,y)``
of the object ``f(x,y)``
from a measured sinogram,
i.e.,
from (usually noisy) samples of ``p(r,\phi)``.


## FBP

A simple image reconstruction method
is called the
"filtered back-projection" (FBP) approach.
It works by filtering each row of the sinogram
with a filter,
called the
[ramp filter](https://en.wikipedia.org/wiki/Radon_transform#Radon_inversion_formula),
whose frequency response
is roughly ``|\nu|``,
where ``\nu`` is the spatial frequency variable
(units cycles / m),
followed by a
[back-projection](https://en.wikipedia.org/wiki/Radon_transform#Dual_transform)
step that is the
[adjoint](https://en.wikipedia.org/wiki/Hermitian_adjoint)
of the Radon transform.


### FBP example

Here is an illustration
of using the `fbp` method
in this package
to perform image reconstruction
from the preceding sinogram.

=#

image = fbp(sino; dr)
(nx,ny) = size(image)
dx = dr # default
x = ((-(nx-1)/2):((nx-1)/2)) * dr
y = x
jim(x, y, image, "FBP image",
    xtick=[-10, 0, 2, 5, 8, 10],
    ytick=[-10, 0, -2, 1, 4, 10],
)


#=

The FBP reconstructed image
looks pretty similar
to a disk of radius ``3``
centered at ``(5,1)``
as expected.
However,
there are quite a few ripples;
these are
[aliasing](https://en.wikipedia.org/wiki/Aliasing)
artifacts
due to the finite sampling
in ``r`` and ``\phi``.

This is example is what is called
2D parallel-beam tomography,
because for the angles ``\phi`` are equally spaced
and for each angle
the radial samples ``r`` are also equally spaced.
This package includes
FBP reconstruction methods
for several other sinogram geometries,
including the well-known fan beam geometries
and the specialized Mojette geometry.


## Noise effects on FBP

Simulating the effects of measurement noise in sinogram
leads to even worse FBP results.

First note that a practical imaging system
has a finite field of view (FOV):
=#

rmax = maximum(r)
fovmask = @. sqrt(abs2(x) + abs2(y)') ≤ rmax
jim(x, y, fovmask, "FOV mask")


# Add noise to the original sinogram:
noisy_sinogram = sino + 0.1 * randn(size(sino))
jimsino(noisy_sinogram, "Noisy sinogram")

# Apply FBP to the noisy sinogram:
noisy_fbp_image = fbp(noisy_sinogram; dr)
noisy_fbp_image .*= fovmask # apply FOV mask
jim(x, y, noisy_fbp_image, "Noisy FBP image"; clim=(0,1))


#=
The methods in this package
(WIP)
and related methods in the
[JuliaImageRecon](https://github.com/JuliaImageRecon)
suite
are designed
to provide better reconstructions
than the simple FBP method.
In particular,
model-based image reconstruction (MBIR) methods
and methods that use suitable machine-learning approaches
can improve image quality significantly.
See this
[2020 survey paper](http://doi.org/10.1109/JPROC.2019.2936204).
=#


# ### Reproducibility

# This page was generated with the following version of Julia:

io = IOBuffer()
versioninfo(io)
split(String(take!(io)), '\n')


# And with the following package versions

import Pkg; Pkg.status()
