# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LinearTransiogram(; ranges, proportions)

A linear transiogram with `ranges`, and `proportions`.

    LinearTransiogram(ball; proportions)

Alternatively, use a custom metric `ball`.
"""
struct LinearTransiogram{B,P} <: Transiogram
  ball::B
  prop::P
end

LinearTransiogram(ball; proportions=(0.5, 0.5)) = LinearTransiogram(ball, proportions)

LinearTransiogram(; ranges=(1.0u"m", 1.0u"m"), proportions=(0.5, 0.5)) =
  LinearTransiogram(MetricBall(ranges), proportions)

constructor(::LinearTransiogram) = LinearTransiogram

function (t::LinearTransiogram)(h)
  r = radius(t.ball)
  p = t.prop
  L = length(p)
  h′, r′ = unitless(h, r)
  v = h′ / r′
  T = typeof(p[1] * v)
  ϵ = T(1e-6) # add small eps for numerical stability
  SMatrix{L,L}(
    i == j ? (h′ < r′) * T(1 - (1 - p[j]) * v) + (h′ ≥ r′) * T(p[j]) - ϵ * (L - 1) :
    (h′ < r′) * T(p[j] * v) + (h′ ≥ r′) * T(p[j]) + ϵ for i in 1:L, j in 1:L
  )
end
