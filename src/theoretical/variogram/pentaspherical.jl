# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PentasphericalVariogram(; range, sill, nugget)
    PentasphericalVariogram(ball; sill, nugget)

A pentaspherical variogram with `range`, `sill` and `nugget`.
Optionally, use a custom metric `ball`.
"""
struct PentasphericalVariogram{V,B} <: Variogram
  sill::V
  nugget::V
  ball::B
  PentasphericalVariogram(sill::V, nugget::V, ball::B) where {V,B} = new{float(V),B}(sill, nugget, ball)
end

PentasphericalVariogram(ball; sill=1.0, nugget=zero(typeof(sill))) = PentasphericalVariogram(sill, nugget, ball)

PentasphericalVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill))) =
  PentasphericalVariogram(sill, nugget, MetricBall(range))

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
