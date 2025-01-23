# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SphericalTransiogram(; ranges, proportions)

A spherical transiogram with `ranges`, and `proportions`.

    SphericalTransiogram(ball; proportions)

Alternatively, use a custom metric `ball`.
"""
struct SphericalTransiogram{B,P} <: Transiogram
  ball::B
  prop::P
end

SphericalTransiogram(ball; proportions=(0.5, 0.5)) = SphericalTransiogram(ball, proportions)

SphericalTransiogram(; ranges=(1.0u"m", 1.0u"m"), proportions=(0.5, 0.5)) =
  SphericalTransiogram(MetricBall(ranges), proportions)

constructor(::SphericalTransiogram) = SphericalTransiogram

function (t::SphericalTransiogram)(h)
  r = radius(t.ball)
  p = t.prop
  L = length(p)
  h′, r′ = unitless(h, r)
  v = 3(h′ / r′) / 2 - (h′ / r′)^3 / 2
  T = typeof(p[1] * v)
  ϵ = T(1e-6) # add small eps for numerical stability
  SMatrix{L,L}(
    i == j ? (h′ < r′) * T(1 - (1 - p[j]) * v) + (h′ ≥ r′) * T(p[j]) - ϵ * (L - 1) :
    (h′ < r′) * T(p[j] * v) + (h′ ≥ r′) * T(p[j]) + ϵ for i in 1:L, j in 1:L
  )
end
