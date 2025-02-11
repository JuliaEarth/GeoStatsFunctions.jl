# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    NuggetEffect(; nugget)

A pure nugget effect variogram with `nugget`.
"""
struct NuggetEffect{V} <: Variogram
  nugget::V
  NuggetEffect(nugget::V) where {V} = new{float(V)}(nugget)
end

NuggetEffect(; nugget=1.0) = NuggetEffect(nugget)

constructor(::NuggetEffect) = NuggetEffect

isstationary(::Type{<:NuggetEffect}) = true

isisotropic(::NuggetEffect) = true

sill(γ::NuggetEffect) = γ.nugget

metricball(::NuggetEffect) = MetricBall(0u"m")

Base.range(::NuggetEffect) = 0u"m"

scale(γ::NuggetEffect, ::Real) = γ

(γ::NuggetEffect)(h) = (h > zero(h)) * γ.nugget

(γ::NuggetEffect)(u::Point, v::Point) = ifelse(u == v, zero(γ.nugget), γ.nugget)
