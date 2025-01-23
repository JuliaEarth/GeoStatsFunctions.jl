# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SineHoleVariogram(; range, sill, nugget)

A sine hole variogram with `range`, `sill` and `nugget`.

    SineHoleVariogram(ball; sill, nugget)

Alternatively, use a custom metric `ball`.
"""
struct SineHoleVariogram{B,V} <: Variogram
  ball::B
  sill::V
  nugget::V
  SineHoleVariogram(ball::B, sill::V, nugget::V) where {B,V} = new{B,float(V)}(ball, sill, nugget)
end

SineHoleVariogram(ball; sill=1.0, nugget=zero(typeof(sill))) = SineHoleVariogram(ball, sill, nugget)

SineHoleVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill))) = SineHoleVariogram(MetricBall(range), sill, nugget)

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
