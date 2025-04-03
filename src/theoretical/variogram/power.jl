# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PowerVariogram(; length, scaling, exponent, nugget)

A power variogram with base `length` in length units,
and `scaling`, `exponent` and `nugget` parameters.

The base `length` parameter serves to scale the lag `h`
in the power variogram formula, i.e. `h -> h / length`.
"""
struct PowerVariogram{ℒ<:Len,V,E} <: Variogram
  length::ℒ
  scaling::V
  nugget::V
  exponent::E
  PowerVariogram(length::ℒ, scaling::V, nugget::V, exponent::E) where {ℒ,V,E} =
    new{float(ℒ),float(V),float(E)}(length, scaling, nugget, exponent)
end

PowerVariogram(; length=1.0, scaling=1.0, nugget=zero(typeof(scaling)), exponent=1.0) =
  PowerVariogram(aslen(length), scaling, nugget, exponent)

constructor(::PowerVariogram) = PowerVariogram

isstationary(::Type{<:PowerVariogram}) = false

isisotropic(::PowerVariogram) = true

metricball(::PowerVariogram) = MetricBall(Inf * u"m")

Base.range(::PowerVariogram) = Inf * u"m"

scale(γ::PowerVariogram, s::Real) = PowerVariogram(s * γ.length, γ.scaling, γ.nugget, γ.exponent)

function (γ::PowerVariogram)(h)
  l = γ.length
  s = γ.scaling
  a = γ.exponent
  n = γ.nugget
  h′, l′ = unitless(h, l)
  s * (h′ / l′)^a + (h′ > 0) * n
end

function (γ::PowerVariogram)(u::Point, v::Point)
  h = evaluate(Euclidean(), u, v)
  γ(h)
end
