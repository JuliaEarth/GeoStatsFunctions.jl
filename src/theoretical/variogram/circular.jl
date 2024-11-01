"""
    CircularVariogram(; range, sill, nugget)
    CircularVariogram(ball; sill, nugget)

A circular variogram with `range`, `sill` and `nugget`.
Optionally, use a custom metric `ball`.
"""
struct CircularVariogram{V,B} <: Variogram
  sill::V
  nugget::V
  ball::B
  CircularVariogram(sill::V, nugget::V, ball::B) where {V,B} = new{float(V),B}(sill, nugget, ball)
end

CircularVariogram(ball; sill=1.0, nugget=zero(typeof(sill))) = CircularVariogram(sill, nugget, ball)

CircularVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill))) = CircularVariogram(sill, nugget, MetricBall(range))

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
