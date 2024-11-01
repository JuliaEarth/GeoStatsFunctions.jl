# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MaternVariogram(; range, sill, nugget, order)
    MaternVariogram(ball; sill, nugget, order)

A Matérn variogram with `range`, `sill` and `nugget`.
Optionally, use a custom metric `ball` and `order` of
the Bessel function.
"""
struct MaternVariogram{V,O,B} <: Variogram
  sill::V
  nugget::V
  order::O
  ball::B
  MaternVariogram(sill::V, nugget::V, order::O, ball::B) where {V,O,B} =
    new{float(V),float(O),B}(sill, nugget, order, ball)
end

MaternVariogram(ball; sill=1.0, nugget=zero(typeof(sill)), order=1.0) = MaternVariogram(sill, nugget, order, ball)

MaternVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill)), order=1.0) =
  MaternVariogram(sill, nugget, order, MetricBall(range))

constructor(::MaternVariogram) = MaternVariogram

isstationary(::Type{<:MaternVariogram}) = true

function (γ::MaternVariogram)(h)
  r = radius(γ.ball)
  s = γ.sill
  n = γ.nugget
  ν = γ.order
  h′, r′ = unitless(h, r)

  # shift lag by machine precision to
  # avoid explosion at the origin
  h′ += eps(typeof(h′))

  δ = √(2ν) * 3(h′ / r′)
  Β = besselk(ν, δ)
  Γ = gamma(ν)

  (s - n) * (1 - 2^(1 - ν) / Γ * δ^ν * Β) + (h′ > 0) * n
end
