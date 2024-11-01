# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GaussianTransiogram(; ranges, proportions)
    GaussianTransiogram(ball; proportions)

A Gaussian transiogram with `ranges`, and `proportions`.
Optionally, use a custom metric `ball`.
"""
struct GaussianTransiogram{P,B} <: Transiogram
  prop::P
  ball::B
end

GaussianTransiogram(ball; proportions=(0.5, 0.5)) = GaussianTransiogram(proportions, ball)

GaussianTransiogram(; ranges=(1.0u"m", 1.0u"m"), proportions=(0.5, 0.5)) =
  GaussianTransiogram(proportions, MetricBall(ranges))

constructor(::GaussianTransiogram) = GaussianTransiogram

function (t::GaussianTransiogram)(h)
  r = radius(t.ball)
  p = t.prop
  L = length(p)
  h′, r′ = unitless(h, r)
  v = 1 - exp(-3(h′ / r′)^2)
  T = typeof(p[1] * v)
  SMatrix{L,L}(i == j ? T(1 - (1 - p[j]) * v) : T((p[j] * v)) for i in 1:L, j in 1:L)
end
