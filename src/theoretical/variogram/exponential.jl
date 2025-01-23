# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExponentialVariogram(; range, sill, nugget)

A exponential variogram with `range`, `sill` and `nugget`.

    ExponentialVariogram(ball; sill, nugget)

Alternatively, use a custom metric `ball`.
"""
struct ExponentialVariogram{B,V} <: Variogram
  ball::B
  sill::V
  nugget::V
  ExponentialVariogram(ball::B, sill::V, nugget::V) where {B,V} = new{B,float(V)}(ball, sill, nugget)
end

ExponentialVariogram(ball; sill=1.0, nugget=zero(typeof(sill))) = ExponentialVariogram(ball, sill, nugget)

ExponentialVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill))) =
  ExponentialVariogram(MetricBall(range), sill, nugget)

constructor(::ExponentialVariogram) = ExponentialVariogram

isstationary(::Type{<:ExponentialVariogram}) = true

function (γ::ExponentialVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  h′, r′ = unitless(h, r)
  (s - n) * (1 - exp(-3(h′ / r′))) + (h′ > 0) * n
end
