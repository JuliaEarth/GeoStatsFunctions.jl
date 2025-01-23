# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GaussianVariogram(; range, sill, nugget)

A Gaussian variogram with `range` in length units,
and `sill` and `nugget` contributions.

    GaussianVariogram(; ranges, rotation, sill, nugget)

Alternatively, use multiple `ranges` and `rotation` matrix
to construct an anisotropic model.

    GaussianVariogram(ball; sill, nugget)

Alternatively, use a custom metric `ball`.

## Examples

```julia
# isotropic model
GaussianVariogram(range=2.0m)

# anisotropic model
GaussianVariogram(ranges=(1.0m, 2.0m))
```
"""
struct GaussianVariogram{B,V} <: Variogram
  ball::B
  sill::V
  nugget::V
  GaussianVariogram(ball::B, sill::V, nugget::V) where {B,V} = new{B,float(V)}(ball, sill, nugget)
end

GaussianVariogram(ball; sill=1.0, nugget=zero(sill)) = GaussianVariogram(ball, sill, nugget)

GaussianVariogram(; range=1.0, ranges=nothing, rotation=I, sill=1.0, nugget=zero(sill)) =
  GaussianVariogram(rangeball(range, ranges, rotation), sill, nugget)

constructor(::GaussianVariogram) = GaussianVariogram

isstationary(::Type{<:GaussianVariogram}) = true

function (γ::GaussianVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget + typeof(s)(1e-6) # add small eps for numerical stability
  h′, r′ = unitless(h, r)
  (s - n) * (1 - exp(-3(h′ / r′)^2)) + (h′ > 0) * n
end
