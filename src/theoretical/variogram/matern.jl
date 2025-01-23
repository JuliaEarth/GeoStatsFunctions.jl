# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MaternVariogram(; range, sill, nugget, order)

A Matérn variogram with `range`, `sill`, `nugget`
and `order` of the Bessel function.

    MaternVariogram(ball; sill, nugget, order)

Alternatively, use a custom metric `ball`.
"""
struct MaternVariogram{B,V,O} <: Variogram
  ball::B
  sill::V
  nugget::V
  order::O
  MaternVariogram(ball::B, sill::V, nugget::V, order::O) where {B,V,O} =
    new{B,float(V),float(O)}(ball, sill, nugget, order)
end

MaternVariogram(ball; sill=1.0, nugget=zero(typeof(sill)), order=1.0) = MaternVariogram(ball, sill, nugget, order)

MaternVariogram(; range=1.0, sill=1.0, nugget=zero(typeof(sill)), order=1.0) =
  MaternVariogram(MetricBall(range), sill, nugget, order)

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
