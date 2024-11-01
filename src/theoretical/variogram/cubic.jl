# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CubicVariogram(; range, sill, nugget)
    CubicVariogram(ball; sill, nugget)

A cubic variogram with `range`, `sill` and `nugget`.
Optionally, use a custom metric `ball`.
"""
struct CubicVariogram{V,B} <: Variogram
  sill::V
  nugget::V
  ball::B
  CubicVariogram(sill::V, nugget::V, ball::B) where {V,B} = new{float(V),B}(sill, nugget, ball)
end

CubicVariogram(ball; sill=1.0, nugget=zero(typeof(sill))) = CubicVariogram(sill, nugget, ball)

CubicVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill))) = CubicVariogram(sill, nugget, MetricBall(range))

constructor(::CubicVariogram) = CubicVariogram

isstationary(::Type{<:CubicVariogram}) = true

function (γ::CubicVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  h′, r′ = unitless(h, r)

  # constants
  T = typeof(h′)
  c1 = T(35) / T(4)
  c2 = T(7) / T(2)
  c3 = T(3) / T(4)
  s1 = 7 * (h′ / r′)^2 - c1 * (h′ / r′)^3 + c2 * (h′ / r′)^5 - c3 * (h′ / r′)^7
  s2 = T(1)

  (h′ < r′) * (s - n) * s1 + (h′ ≥ r′) * (s - n) * s2 + (h′ > 0) * n
end
