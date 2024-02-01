# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    NuggetEffect(nugget=n)

A pure nugget effect variogram with nugget `n`.
"""
struct NuggetEffect{V} <: Variogram
  nugget::V
end

NuggetEffect(; nugget=1.0) = NuggetEffect(nugget)

(γ::NuggetEffect)(h) = (h > 0) * γ.nugget

(γ::NuggetEffect)(u::Point, v::Point) = ifelse(u == v, zero(γ.nugget), γ.nugget)

Base.range(::NuggetEffect{T}) where {T} = zero(T)

sill(γ::NuggetEffect) = γ.nugget

variotype(::NuggetEffect) = NuggetEffect

isstationary(::Type{<:NuggetEffect}) = true

isisotropic(::NuggetEffect) = true
