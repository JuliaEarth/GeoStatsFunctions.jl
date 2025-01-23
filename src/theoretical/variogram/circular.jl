# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CircularVariogram(; range, sill, nugget)

A circular variogram with `range` in length units,
and `sill` and `nugget` contributions.

    CircularVariogram(; ranges, rotation, sill, nugget)

Alternatively, use multiple `ranges` and `rotation` matrix
to construct an anisotropic model.

    CircularVariogram(ball; sill, nugget)

Alternatively, use a custom metric `ball`.

## Examples

```julia
# isotropic model
CircularVariogram(range=2.0m)

# anisotropic model
CircularVariogram(ranges=(1.0m, 2.0m))
```
"""
struct CircularVariogram{B,V} <: Variogram
  ball::B
  sill::V
  nugget::V
  CircularVariogram(ball::B, sill::V, nugget::V) where {B,V} = new{B,float(V)}(ball, sill, nugget)
end

CircularVariogram(ball; sill=1.0, nugget=zero(sill)) = CircularVariogram(ball, sill, nugget)

CircularVariogram(; range=1.0, ranges=nothing, rotation=I, sill=1.0, nugget=zero(sill)) =
  CircularVariogram(rangeball(range, ranges, rotation), sill, nugget)

constructor(::CircularVariogram) = CircularVariogram

isstationary(::Type{<:CircularVariogram}) = true

function (γ::CircularVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  h′, r′ = unitless(h, r)
  v = h′ ≤ r′ ? 1 - (2 / π) * acos(h′ / r′) + (2h′ / (π * r′)) * sqrt(1 - (h′^2 / r′^2)) : one(h′)
  (s - n) * v + (h′ > zero(h′)) * n
end
