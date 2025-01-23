"""
    CircularVariogram(; range, sill, nugget)

A circular variogram with `range`, `sill` and `nugget`.

    CircularVariogram(ball; sill, nugget)

Alternatively, use a custom metric `ball`.
"""
struct CircularVariogram{B,V} <: Variogram
  ball::B
  sill::V
  nugget::V
  CircularVariogram(ball::B, sill::V, nugget::V) where {B,V} = new{B,float(V)}(ball, sill, nugget)
end

CircularVariogram(ball; sill=1.0, nugget=zero(typeof(sill))) = CircularVariogram(ball, sill, nugget)

CircularVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill))) = CircularVariogram(MetricBall(range), sill, nugget)

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
