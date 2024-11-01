# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PowerVariogram(; scaling, exponent, nugget)

A power variogram with `scaling`, `exponent` and `nugget`.
"""
struct PowerVariogram{V,E} <: Variogram
  scaling::V
  nugget::V
  exponent::E
  PowerVariogram(scaling::V, nugget::V, exponent::E) where {V,E} = new{float(V),float(E)}(scaling, nugget, exponent)
end

PowerVariogram(; scaling=1.0, nugget=zero(typeof(scaling)), exponent=1.0) = PowerVariogram(scaling, nugget, exponent)

constructor(::PowerVariogram) = PowerVariogram

isstationary(::Type{<:PowerVariogram}) = false

isisotropic(::PowerVariogram) = true

function (γ::PowerVariogram)(h)
  s = γ.scaling
  a = γ.exponent
  n = γ.nugget
  h′ = ustrip(h)
  s * h′^a + (h′ > 0) * n
end

function (γ::PowerVariogram)(u::Point, v::Point)
  d = Euclidean()
  h = evaluate(d, u, v)
  γ(h)
end
