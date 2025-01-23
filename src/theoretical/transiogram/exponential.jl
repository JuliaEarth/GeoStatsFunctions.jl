# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExponentialTransiogram(; ranges, proportions)

A exponential transiogram with `ranges`, and `proportions`.

    ExponentialTransiogram(ball; proportions)

Alternatively, use a custom metric `ball`.
"""
struct ExponentialTransiogram{B,P} <: Transiogram
  ball::B
  prop::P
end

ExponentialTransiogram(ball; proportions=(0.5, 0.5)) = ExponentialTransiogram(ball, proportions)

ExponentialTransiogram(; ranges=(1.0u"m", 1.0u"m"), proportions=(0.5, 0.5)) =
  ExponentialTransiogram(MetricBall(ranges), proportions)

constructor(::ExponentialTransiogram) = ExponentialTransiogram

function (t::ExponentialTransiogram)(h)
  r = radius(t.ball)
  p = t.prop
  L = length(p)
  h′, r′ = unitless(h, r)
  v = 1 - exp(-3(h′ / r′))
  T = typeof(p[1] * v)
  ϵ = T(1e-6) # add small eps for numerical stability
  SMatrix{L,L}(i == j ? T(1 - (1 - p[j]) * v) - ϵ * (L - 1) : T(p[j] * v) + ϵ for i in 1:L, j in 1:L)
end
