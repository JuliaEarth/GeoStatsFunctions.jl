# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExponentialVariogram(range=r, sill=s, nugget=n)
    ExponentialVariogram(ball; sill=s, nugget=n)

An exponential variogram with range `r`, sill `s` and nugget `n`.
Optionally, use a custom metric `ball`.
"""
struct ExponentialVariogram{V,B} <: Variogram
  sill::V
  nugget::V
  ball::B
  ExponentialVariogram(sill::V, nugget::V, ball::B) where {V,B} = new{float(V),B}(sill, nugget, ball)
end

ExponentialVariogram(ball; sill=1.0, nugget=zero(typeof(sill))) = ExponentialVariogram(sill, nugget, ball)

ExponentialVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill))) =
  ExponentialVariogram(sill, nugget, MetricBall(range))

constructor(::ExponentialVariogram) = ExponentialVariogram

isstationary(::Type{<:ExponentialVariogram}) = true

function (γ::ExponentialVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  h′, r′ = unitless(h, r)
  (s - n) * (1 - exp(-3(h′ / r′))) + (h′ > 0) * n
end
