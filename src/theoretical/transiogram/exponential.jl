# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExponentialTransiogram(; range, proportions)

An exponential transiogram with `range` in length units,
and categorical `proportions`.

    ExponentialTransiogram(; ranges, rotation, proportions)

Alternatively, use multiple `ranges` and `rotation` matrix
to construct an anisotropic model.

    ExponentialTransiogram(ball; proportions)

Alternatively, use a custom metric `ball`.
"""
struct ExponentialTransiogram{B,P} <: Transiogram
  ball::B
  proportions::P
end

ExponentialTransiogram(ball; proportions=(0.5, 0.5)) = ExponentialTransiogram(ball, proportions)

ExponentialTransiogram(; range=1.0, ranges=nothing, rotation=I, proportions=(0.5, 0.5)) =
  ExponentialTransiogram(rangeball(range, ranges, rotation), proportions)

constructor(::ExponentialTransiogram) = ExponentialTransiogram

function (t::ExponentialTransiogram)(h)
  r = radius(t.ball)
  p = t.proportions
  L = length(p)
  h′, r′ = unitless(h, r)
  v = 1 - exp(-3(h′ / r′))
  T = typeof(p[1] * v)
  SMatrix{L,L}(i == j ? T(1 - (1 - p[j]) * v) : T(p[j] * v) for i in 1:L, j in 1:L)
end
