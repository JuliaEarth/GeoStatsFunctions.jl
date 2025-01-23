# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SineHoleVariogram(; range, sill, nugget)

A sinehole variogram with `range` in length units,
and `sill` and `nugget` contributions.

    SineHoleVariogram(; ranges, rotation, sill, nugget)

Alternatively, use multiple `ranges` and `rotation` matrix
to construct an anisotropic model.

    SineHoleVariogram(ball; sill, nugget)

Alternatively, use a custom metric `ball`.

## Examples

```julia
# isotropic model
SineHoleVariogram(range=2.0m)

# anisotropic model
SineHoleVariogram(ranges=(1.0m, 2.0m))
```
"""
struct SineHoleVariogram{B,V} <: Variogram
  ball::B
  sill::V
  nugget::V
  SineHoleVariogram(ball::B, sill::V, nugget::V) where {B,V} = new{B,float(V)}(ball, sill, nugget)
end

SineHoleVariogram(ball; sill=1.0, nugget=zero(sill)) = SineHoleVariogram(ball, sill, nugget)

SineHoleVariogram(; range=1.0, ranges=nothing, rotation=I, sill=1.0, nugget=zero(sill)) =
  SineHoleVariogram(rangeball(range, ranges, rotation), sill, nugget)

constructor(::SineHoleVariogram) = SineHoleVariogram

isstationary(::Type{<:SineHoleVariogram}) = true

function (γ::SineHoleVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  h′, r′ = unitless(h, r)

  # shift lag by machine precision to
  # avoid explosion at the origin
  h′ += eps(typeof(h′))
  c = oftype(h′, π)

  (s - n) * (1 - sin(c * h′ / r′) / (c * h′ / r′)) + (h′ > 0) * n
end
