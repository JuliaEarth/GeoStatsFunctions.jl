# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const Len{T} = Quantity{T,u"𝐋"}
const InvLen{T} = Quantity{T,u"𝐋^-1"}

aslen(x::Len) = x
aslen(x::Number) = x * u"m"
aslen(::Quantity) = throw(ArgumentError("invalid length unit"))

asinvlen(x::InvLen) = x
asinvlen(x::Number) = x * u"m^-1"
asinvlen(::Quantity) = throw(ArgumentError("invalid length unit"))

unitless(a::Number, b::Quantity) = a, ustrip(b)
function unitless(a::Quantity, b::Quantity)
  u = Unitful.promote_unit(unit(a), unit(b))
  ustrip(u, a), ustrip(u, b)
end

isinvalid(v) = ismissing(v) || (v isa Number && !isfinite(v))

defaultmaxlag(data) = minside(boundingbox(domain(data))) / 2

defaultvariables(n) = n == 1 ? [:field] : [Symbol(:field, i) for i in 1:n]

rangeball(range, ::Nothing, rotation) = MetricBall(range)
rangeball(range, ranges, rotation) = MetricBall(ranges, rotation)

function minside(box)
  s = uniquesides(box)
  minimum(filter(>(zero(eltype(s))), s))
end

uniquesides(box::Box{<:𝔼}) = sides(box)

function uniquesides(box::Box{<:🌐})
  r = vertices(boundary(box))
  s1 = length(Segment(r[1], r[2]))
  s2 = length(Segment(r[2], r[3]))
  s3 = length(Segment(r[3], r[4]))
  s4 = length(Segment(r[4], r[1]))
  (s1, s2, s3, s4)
end

function printlnvec(io, vec, n)
  printvec(io, vec, n)
  println(io)
end

function printvec(io, vec::AbstractArray, n)
  print(io, "[")
  if length(vec) > 2n
    k = n - 1
    join(io, vec[begin:(begin + k)], ", ")
    print(io, ", ..., ")
    join(io, vec[(end - k):end], ", ")
  else
    join(io, vec, ", ")
  end
  print(io, "]")
end

function printvec(io, vec::AbstractArray{<:AbstractArray}, n)
  len = length(vec)
  println(io)
  if len > 2n
    for i in 1:n
      print(io, "│  ├─ ")
      printlnvec(io, vec[i], n)
    end
    println(io, "│  ⋮")
    for i in (len - n + 1):(len - 1)
      print(io, "│  ├─ ")
      printlnvec(io, vec[i], n)
    end
  else
    for i in 1:(len - 1)
      print(io, "│  ├─ ")
      printlnvec(io, vec[i], n)
    end
  end
  print(io, "│  └─ ")
  printvec(io, vec[len], n)
end
