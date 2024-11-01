# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LinearTransiogram(; ranges, proportions)
    LinearTransiogram(ball; proportions)

A linear transiogram with `ranges`, and `proportions`.
Optionally, use a custom metric `ball`.
"""
struct LinearTransiogram{P,B} <: Transiogram
  prop::P
  ball::B
end

LinearTransiogram(ball; proportions=(0.5, 0.5)) = LinearTransiogram(proportions, ball)

LinearTransiogram(; ranges=(1.0u"m", 1.0u"m"), proportions=(0.5, 0.5)) =
  LinearTransiogram(proportions, MetricBall(ranges))

constructor(::LinearTransiogram) = LinearTransiogram

function (t::LinearTransiogram)(h)
  r = radius(t.ball)
  p = t.prop
  L = length(p)
  h′, r′ = unitless(h, r)
  SMatrix{L,L}(
    i == j ? (h′ < p[j]) * (1 - (1 - p[j]) * (h′ / r′)) + (h′ ≥ p[j]) * p[j] :
    (h′ < p[j]) * (p[j] * (h′ / r′)) + (h′ ≥ p[j]) * p[j] for i in 1:L, j in 1:L
  )
end
