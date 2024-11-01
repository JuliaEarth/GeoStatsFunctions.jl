# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExponentialTransiogram(; ranges, proportions)
    ExponentialTransiogram(ball; proportions)

A exponential transiogram with `ranges`, and `proportions`.
Optionally, use a custom metric `ball`.
"""
struct ExponentialTransiogram{P,B} <: Transiogram
  prop::P
  ball::B
end

ExponentialTransiogram(ball; proportions=(0.5, 0.5)) = ExponentialTransiogram(proportions, ball)

ExponentialTransiogram(; ranges=(1.0u"m", 1.0u"m"), proportions=(0.5, 0.5)) =
  ExponentialTransiogram(proportions, MetricBall(ranges))

constructor(::ExponentialTransiogram) = ExponentialTransiogram

function (t::ExponentialTransiogram)(h)
  r = radius(t.ball)
  p = t.prop
  L = length(p)
  h′, r′ = unitless(h, r)
  v = 1 - exp(-3(h′ / r′))
  T = typeof(p[1] * v)
  SMatrix{L,L}(i == j ? T(1 - (1 - p[j]) * v) : T((p[j] * v)) for i in 1:L, j in 1:L)
end
