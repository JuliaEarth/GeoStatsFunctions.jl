# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SineHoleVariogram(; range, sill, nugget)
    SineHoleVariogram(ball; sill, nugget)

A sine hole variogram with `range`, `sill` and `nugget`.
Optionally, use a custom metric `ball`.
"""
struct SineHoleVariogram{V,B} <: Variogram
  sill::V
  nugget::V
  ball::B
  SineHoleVariogram(sill::V, nugget::V, ball::B) where {V,B} = new{float(V),B}(sill, nugget, ball)
end

SineHoleVariogram(ball; sill=1.0, nugget=zero(typeof(sill))) = SineHoleVariogram(sill, nugget, ball)

SineHoleVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill))) = SineHoleVariogram(sill, nugget, MetricBall(range))

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
