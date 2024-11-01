# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SphericalTransiogram(; ranges, proportions)
    SphericalTransiogram(ball; proportions)

A spherical transiogram with `ranges`, and `proportions`.
Optionally, use a custom metric `ball`.
"""
struct SphericalTransiogram{P,B} <: Transiogram
  prop::P
  ball::B
end

SphericalTransiogram(ball; proportions=(0.5, 0.5)) = SphericalTransiogram(proportions, ball)

SphericalTransiogram(; ranges=(1.0u"m", 1.0u"m"), proportions=(0.5, 0.5)) =
  SphericalTransiogram(proportions, MetricBall(ranges))

constructor(::SphericalTransiogram) = SphericalTransiogram

function (t::SphericalTransiogram)(h)
  r = radius(t.ball)
  p = t.prop
  L = length(p)
  h′, r′ = unitless(h, r)
  v = 3(h′ / r′) / 2 - (h′ / r′)^3 / 2
  T = typeof(p[1] * v)
  SMatrix{L,L}(
    i == j ? (h′ < p[j]) * T(1 - (1 - p[j]) * v) + (h′ ≥ p[j]) * T(p[j]) :
    (h′ < p[j]) * T((p[j] * v)) + (h′ ≥ p[j]) * T(p[j]) for i in 1:L, j in 1:L
  )
end
