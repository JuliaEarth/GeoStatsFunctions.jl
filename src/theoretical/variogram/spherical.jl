# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SphericalVariogram(; range, sill, nugget)

A spherical variogram with `range` in length units,
and `sill` and `nugget` contributions.

    SphericalVariogram(; ranges, rotation, sill, nugget)

Alternatively, use multiple `ranges` and `rotation` matrix
to construct an anisotropic model.

    SphericalVariogram(ball; sill, nugget)

Alternatively, use a custom metric `ball`.

## Examples

```julia
# isotropic model
SphericalVariogram(range=2.0m)

# anisotropic model
SphericalVariogram(ranges=(1.0m, 2.0m))
```
"""
struct SphericalVariogram{B,V} <: Variogram
  ball::B
  sill::V
  nugget::V
  SphericalVariogram(ball::B, sill::V, nugget::V) where {B,V} = new{B,float(V)}(ball, sill, nugget)
end

SphericalVariogram(ball; sill=1.0, nugget=zero(sill)) = SphericalVariogram(ball, sill, nugget)

SphericalVariogram(; range=1.0, ranges=nothing, rotation=I, sill=1.0, nugget=zero(sill)) =
  SphericalVariogram(rangeball(range, ranges, rotation), sill, nugget)

constructor(::SphericalVariogram) = SphericalVariogram

isstationary(::Type{<:SphericalVariogram}) = true

function (γ::SphericalVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  h′, r′ = unitless(h, r)

  # constants
  T = typeof(h′)
  c1 = T(3) / T(2)
  c2 = T(1) / T(2)
  s1 = c1 * (h′ / r′) - c2 * (h′ / r′)^3
  s2 = T(1)

  (h′ < r′) * (s - n) * s1 + (h′ ≥ r′) * (s - n) * s2 + (h′ > 0) * n
end
