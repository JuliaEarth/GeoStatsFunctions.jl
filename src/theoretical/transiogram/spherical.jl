# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SphericalTransiogram(; range, proportions)

A spherical transiogram with `range` in length units,
and categorical `proportions`.

    SphericalTransiogram(; ranges, rotation, proportions)

Alternatively, use multiple `ranges` and `rotation` matrix
to construct an anisotropic model.

    SphericalTransiogram(ball; proportions)

Alternatively, use a custom metric `ball`.
"""
struct SphericalTransiogram{B,P<:NTuple} <: Transiogram
  ball::B
  proportions::P
end

SphericalTransiogram(ball; proportions=(0.5, 0.5)) = SphericalTransiogram(ball, proportions)

SphericalTransiogram(; range=1.0, ranges=nothing, rotation=I, proportions=(0.5, 0.5)) =
  SphericalTransiogram(rangeball(range, ranges, rotation), proportions)

constructor(::SphericalTransiogram) = SphericalTransiogram

function (t::SphericalTransiogram)(h)
  r = radius(t.ball)
  p = t.proportions
  L = length(p)
  h′, r′ = unitless(h, r)
  v = 3(h′ / r′) / 2 - (h′ / r′)^3 / 2
  T = typeof(p[1] * v)
  SMatrix{L,L}(
    i == j ? (h′ < r′) * T(1 - (1 - p[j]) * v) + (h′ ≥ r′) * T(p[j]) : (h′ < r′) * T(p[j] * v) + (h′ ≥ r′) * T(p[j]) for
    i in 1:L, j in 1:L
  )
end
