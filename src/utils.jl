# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const Len{T} = Quantity{T,u"𝐋"}
const InvLen{T} = Quantity{T,u"𝐋^-1"}

asinvlen(x::Number) = x * u"m^-1"
asinvlen(x::InvLen) = x
asinvlen(::Quantity) = throw(ArgumentError("invalid units, please check the documentation"))

addunit(x::Number, u) = x * u
addunit(::Quantity, _) = throw(ArgumentError("invalid units, please check the documentation"))

unitless(a::Number, b::Quantity) = a, ustrip(b)
function unitless(a::Quantity, b::Quantity)
  u = Unitful.promote_unit(unit(a), unit(b))
  ustrip(u, a), ustrip(u, b)
end

defaultmaxlag(data) = _minside(boundingbox(domain(data))) / 2

function _minside(box)
  s = _sides(box)
  minimum(filter(>(zero(eltype(s))), s))
end

_sides(box::Box{<:𝔼}) = sides(box)

function _sides(box::Box{<:🌐})
  r = vertices(boundary(box))
  s1 = length(Segment(r[1], r[2]))
  s2 = length(Segment(r[2], r[3]))
  s3 = length(Segment(r[3], r[4]))
  s4 = length(Segment(r[4], r[1]))
  (s1, s2, s3, s4)
end
