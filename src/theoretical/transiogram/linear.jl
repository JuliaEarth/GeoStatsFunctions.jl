# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LinearTransiogram(; range, proportions)

A linear transiogram with `range` in length units,
and categorical `proportions`.

    LinearTransiogram(; ranges, rotation, proportions)

Alternatively, use multiple `ranges` and `rotation` matrix
to construct an anisotropic model.

    LinearTransiogram(ball; proportions)

Alternatively, use a custom metric `ball`.
"""
struct LinearTransiogram{B,P<:NTuple} <: Transiogram
  ball::B
  proportions::P
end

LinearTransiogram(ball; proportions=(0.5, 0.5)) = LinearTransiogram(ball, proportions)

LinearTransiogram(; range=1.0, ranges=nothing, rotation=I, proportions=(0.5, 0.5)) =
  LinearTransiogram(rangeball(range, ranges, rotation), proportions)

constructor(::LinearTransiogram) = LinearTransiogram

function (t::LinearTransiogram)(h)
  r = radius(t.ball)
  p = t.proportions
  L = length(p)
  h′, r′ = unitless(h, r)
  v = h′ / r′
  T = typeof(p[1] * v)
  SMatrix{L,L}(
    i == j ? (h′ < r′) * T(1 - (1 - p[j]) * v) + (h′ ≥ r′) * T(p[j]) : (h′ < r′) * T(p[j] * v) + (h′ ≥ r′) * T(p[j]) for
    i in 1:L, j in 1:L
  )
end
