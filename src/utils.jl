# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const Len{T} = Quantity{T,u"𝐋"}

addunit(x::Number, u) = x * u
addunit(::Quantity, _) = throw(ArgumentError("invalid units, please check the documentation"))

unitless(a::Number, b::Quantity) = a, ustrip(b)
function unitless(a::Quantity, b::Quantity)
  u = Unitful.promote_unit(unit(a), unit(b))
  ustrip(u, a), ustrip(u, b)
end

function uevaluate(distance, p₁, p₂)
  u₁ = unit(Meshes.lentype(p₁))
  u₂ = unit(Meshes.lentype(p₂))
  u = Unitful.promote_unit(u₁, u₂)
  v₁ = ustrip.(u, to(p₁))
  v₂ = ustrip.(u, to(p₂))
  evaluate(distance, v₁, v₂) * u
end
