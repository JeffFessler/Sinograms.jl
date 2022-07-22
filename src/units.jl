# units.jl

using .Unitful: °, rad

"""
    to_radians(angle::Unitful.Quantity)

Convert Unitful quantity to radians.
"""
to_radians(angle::Unitful.Quantity) = convert(eltype(1.0rad), angle)
