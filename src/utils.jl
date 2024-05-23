# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const Len{T} = Quantity{T,u"𝐋"}

addunit(x::Number, u) = x * u
addunit(::Quantity, _) = throw(ArgumentError("invalid units, please check the documentation"))

unitless(xs::Quantity...) = ustrip.(promote(xs...))

function uevaluate(distance, p₁, p₂)
  u₁ = unit(Meshes.lentype(p₁))
  u₂ = unit(Meshes.lentype(p₂))
  u = Unitful.promote_unit(u₁, u₂)
  v₁ = ustrip.(u, to(p₁))
  v₂ = ustrip.(u, to(p₂))
  evaluate(distance, v₁, v₂) * u
end

"""
    spheredir(θ, φ)

Return the 3D direction given polar angle `θ` and
azimuthal angle `φ` in degrees according to the ISO
convention.
"""
function spheredir(theta, phi)
  θ, φ = deg2rad(theta), deg2rad(phi)
  Vec(sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ))
end
