# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExponentialVariogram(; range, sill, nugget)

An exponential variogram with `range` in length units,
and `sill` and `nugget` contributions.

    ExponentialVariogram(; ranges, rotation, sill, nugget)

Alternatively, use multiple `ranges` and `rotation` matrix
to construct an anisotropic model.

    ExponentialVariogram(ball; sill, nugget)

Alternatively, use a custom metric `ball`.

## Examples

```julia
# isotropic model
ExponentialVariogram(range=2.0m)

# anisotropic model
ExponentialVariogram(ranges=(1.0m, 2.0m))
```
"""
struct ExponentialVariogram{B,V} <: Variogram
  ball::B
  sill::V
  nugget::V
  ExponentialVariogram(ball::B, sill::V, nugget::V) where {B,V} = new{B,float(V)}(ball, sill, nugget)
end

ExponentialVariogram(ball; sill=1.0, nugget=zero(sill)) = ExponentialVariogram(ball, sill, nugget)

ExponentialVariogram(; range=1.0, ranges=nothing, rotation=I, sill=1.0, nugget=zero(sill)) =
  ExponentialVariogram(rangeball(range, ranges, rotation), sill, nugget)

constructor(::ExponentialVariogram) = ExponentialVariogram

isstationary(::Type{<:ExponentialVariogram}) = true

function (γ::ExponentialVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  h′, r′ = unitless(h, r)
  (s - n) * (1 - exp(-3(h′ / r′))) + (h′ > 0) * n
end
