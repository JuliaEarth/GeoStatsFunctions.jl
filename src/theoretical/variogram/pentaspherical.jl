# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PentaSphericalVariogram(; range, sill, nugget)

A pentaspherical variogram with `range` in length units,
and `sill` and `nugget` contributions.

    PentaSphericalVariogram(; ranges, rotation, sill, nugget)

Alternatively, use multiple `ranges` and `rotation` matrix
to construct an anisotropic model.

    PentaSphericalVariogram(ball; sill, nugget)

Alternatively, use a custom metric `ball`.

## Examples

```julia
# isotropic model
PentaSphericalVariogram(range=2.0m)

# anisotropic model
PentaSphericalVariogram(ranges=(1.0m, 2.0m))
```
"""
struct PentaSphericalVariogram{B,V} <: Variogram
  ball::B
  sill::V
  nugget::V
  PentaSphericalVariogram(ball::B, sill::V, nugget::V) where {B,V} = new{B,float(V)}(ball, sill, nugget)
end

PentaSphericalVariogram(ball; sill=1.0, nugget=zero(sill)) = PentaSphericalVariogram(ball, sill, nugget)

PentaSphericalVariogram(; range=1.0, ranges=nothing, rotation=I, sill=1.0, nugget=zero(sill)) =
  PentaSphericalVariogram(rangeball(range, ranges, rotation), sill, nugget)

constructor(::PentaSphericalVariogram) = PentaSphericalVariogram

isstationary(::Type{<:PentaSphericalVariogram}) = true

function (γ::PentaSphericalVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  h′, r′ = unitless(h, r)

  # constants
  T = typeof(h′)
  c1 = T(15) / T(8)
  c2 = T(5) / T(4)
  c3 = T(3) / T(8)
  s1 = c1 * (h′ / r′) - c2 * (h′ / r′)^3 + c3 * (h′ / r′)^5
  s2 = T(1)

  (h′ < r′) * (s - n) * s1 + (h′ ≥ r′) * (s - n) * s2 + (h′ > 0) * n
end
