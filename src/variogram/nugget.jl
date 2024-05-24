# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    NuggetEffect(nugget=n)

A pure nugget effect variogram with nugget `n`.
"""
struct NuggetEffect{V} <: Variogram
  nugget::V
end

NuggetEffect(; nugget=1.0) = NuggetEffect(nugget)

variotype(::NuggetEffect) = NuggetEffect

isstationary(::Type{<:NuggetEffect}) = true

isisotropic(::NuggetEffect) = true

sill(γ::NuggetEffect) = γ.nugget

Base.range(::NuggetEffect) = 0u"m"

scale(γ::NuggetEffect, ::Real) = γ

(γ::NuggetEffect)(h) = (h > zero(h)) * γ.nugget

(γ::NuggetEffect)(u::Point, v::Point) = ifelse(u == v, zero(γ.nugget), γ.nugget)
