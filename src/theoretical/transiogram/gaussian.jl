# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GaussianTransiogram(; range, proportions)

A Gaussian transiogram with `range` in length units,
and categorical `proportions`.

    GaussianTransiogram(; ranges, rotation, proportions)

Alternatively, use multiple `ranges` and `rotation` matrix
to construct an anisotropic model.

    GaussianTransiogram(ball; proportions)

Alternatively, use a custom metric `ball`.
"""
struct GaussianTransiogram{B,P<:NTuple} <: Transiogram
  ball::B
  proportions::P
end

GaussianTransiogram(ball; proportions=(0.5, 0.5)) = GaussianTransiogram(ball, proportions)

GaussianTransiogram(; range=1.0, ranges=nothing, rotation=I, proportions=(0.5, 0.5)) =
  GaussianTransiogram(rangeball(range, ranges, rotation), proportions)

constructor(::GaussianTransiogram) = GaussianTransiogram

function (t::GaussianTransiogram)(h)
  r = radius(t.ball)
  p = t.proportions
  L = length(p)
  h′, r′ = unitless(h, r)
  v = 1 - exp(-3(h′ / r′)^2)
  T = typeof(p[1] * v)
  SMatrix{L,L}(i == j ? T(1 - (1 - p[j]) * v) : T(p[j] * v) for i in 1:L, j in 1:L)
end
