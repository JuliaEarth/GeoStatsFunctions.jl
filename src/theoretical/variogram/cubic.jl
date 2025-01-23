# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CubicVariogram(; range, sill, nugget)

A cubic variogram with `range` in length units,
and `sill` and `nugget` contributions.

    CubicVariogram(; ranges, rotation, sill, nugget)

Alternatively, use multiple `ranges` and `rotation` matrix
to construct an anisotropic model.

    CubicVariogram(ball; sill, nugget)

Alternatively, use a custom metric `ball`.

## Examples

```julia
# isotropic model
CubicVariogram(range=2.0m)

# anisotropic model
CubicVariogram(ranges=(1.0m, 2.0m))
```
"""
struct CubicVariogram{B,V} <: Variogram
  ball::B
  sill::V
  nugget::V
  CubicVariogram(ball::B, sill::V, nugget::V) where {B,V} = new{B,float(V)}(ball, sill, nugget)
end

CubicVariogram(ball; sill=1.0, nugget=zero(sill)) = CubicVariogram(ball, sill, nugget)

CubicVariogram(; range=1.0, ranges=nothing, rotation=I, sill=1.0, nugget=zero(sill)) =
  CubicVariogram(rangeball(range, ranges, rotation), sill, nugget)

constructor(::CubicVariogram) = CubicVariogram

isstationary(::Type{<:CubicVariogram}) = true

function (γ::CubicVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  h′, r′ = unitless(h, r)

  # constants
  T = typeof(h′)
  c1 = T(35) / T(4)
  c2 = T(7) / T(2)
  c3 = T(3) / T(4)
  s1 = 7 * (h′ / r′)^2 - c1 * (h′ / r′)^3 + c2 * (h′ / r′)^5 - c3 * (h′ / r′)^7
  s2 = T(1)

  (h′ < r′) * (s - n) * s1 + (h′ ≥ r′) * (s - n) * s2 + (h′ > 0) * n
end
