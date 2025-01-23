# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GaussianTransiogram(; ranges, proportions)

A Gaussian transiogram with `ranges`, and `proportions`.

    GaussianTransiogram(ball; proportions)

Alternatively, use a custom metric `ball`.
"""
struct GaussianTransiogram{B,P} <: Transiogram
  ball::B
  prop::P
end

GaussianTransiogram(ball; proportions=(0.5, 0.5)) = GaussianTransiogram(ball, proportions)

GaussianTransiogram(; ranges=(1.0u"m", 1.0u"m"), proportions=(0.5, 0.5)) =
  GaussianTransiogram(MetricBall(ranges), proportions)

constructor(::GaussianTransiogram) = GaussianTransiogram

function (t::GaussianTransiogram)(h)
  r = radius(t.ball)
  p = t.prop
  L = length(p)
  h′, r′ = unitless(h, r)
  v = 1 - exp(-3(h′ / r′)^2)
  T = typeof(p[1] * v)
  ϵ = T(1e-6) # add small eps for numerical stability
  SMatrix{L,L}(i == j ? T(1 - (1 - p[j]) * v) - ϵ * (L - 1) : T(p[j] * v) + ϵ for i in 1:L, j in 1:L)
end
