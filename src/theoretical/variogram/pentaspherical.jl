# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PentasphericalVariogram(; range, sill, nugget)

A pentaspherical variogram with `range`, `sill` and `nugget`.

    PentasphericalVariogram(ball; sill, nugget)

Alternatively, use a custom metric `ball`.
"""
struct PentasphericalVariogram{B,V} <: Variogram
  ball::B
  sill::V
  nugget::V
  PentasphericalVariogram(ball::B, sill::V, nugget::V) where {B,V} = new{B,float(V)}(ball, sill, nugget)
end

PentasphericalVariogram(ball; sill=1.0, nugget=zero(typeof(sill))) = PentasphericalVariogram(ball, sill, nugget)

PentasphericalVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill))) =
  PentasphericalVariogram(MetricBall(range), sill, nugget)

constructor(::PentasphericalVariogram) = PentasphericalVariogram

isstationary(::Type{<:PentasphericalVariogram}) = true

function (γ::PentasphericalVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  h′, r′ = unitless(h, r)

  # constants
  T = typeof(h′)
  c1 = T(15) / T(8)
  c2 = T(5) / T(4)
  c3 = T(3) / T(8)
  s1 = c1 * (h′ / r′) - c2 * (h′ / r′)^3 + c3 * (h′ / r′)^5
  s2 = T(1)

  (h′ < r′) * (s - n) * s1 + (h′ ≥ r′) * (s - n) * s2 + (h′ > 0) * n
end
