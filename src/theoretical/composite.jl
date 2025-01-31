# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# helper function to extract raw data
# from uniform scaling objects
raw(a::UniformScaling) = a.λ
raw(a) = a

"""
    CompositeFunction(cs, fs)

A composite geostatistical function `f = c₁f₁ + c₂f₂ + ⋯ + cₙfₙ` with
coefficients `cs = (c₁, c₂, ..., cₙ)` and geostatistical functions
`fs = (f₁, f₂, ..., fₙ)`.
"""
struct CompositeFunction{CS,FS} <: GeoStatsFunction
  cs::CS
  fs::FS
end

isstationary(f::CompositeFunction) = all(isstationary, f.fs)

isisotropic(f::CompositeFunction) = all(isisotropic, f.fs)

isbanded(f::CompositeFunction) = all(isbanded, f.fs)

sill(f::CompositeFunction) = raw(sum(f.cs .* map(sill, f.fs)))

nugget(f::CompositeFunction) = raw(sum(f.cs .* map(nugget, f.fs)))

Base.range(f::CompositeFunction) = maximum(range(fᵢ) for fᵢ in f.fs)

scale(f::CompositeFunction, s::Real) = CompositeFunction(f.cs, map(fᵢ -> scale(fᵢ, s), f.fs))

function structures(f::CompositeFunction)
  ks, fs = f.cs, f.fs

  # total nugget and contributions
  cₒ = raw(sum(@. ks * nugget(fs)))
  cs = @. raw(ks * (sill(fs) - nugget(fs)))

  # discard nugget effect terms
  inds = findall(fᵢ -> !(fᵢ isa NuggetEffect), fs)
  cs, fs = cs[inds], fs[inds]

  # adjust sill and nugget
  fs = map(f -> first(last(structures(f))), fs)

  cₒ, cs, fs
end

(f::CompositeFunction)(h) = raw(sum(f.cs .* map(fᵢ -> fᵢ(h), f.fs)))
(f::CompositeFunction)(u::Point, v::Point) = raw(sum(f.cs .* map(fᵢ -> fᵢ(u, v), f.fs)))

# algebraic structure
*(c, f::GeoStatsFunction) = CompositeFunction((c,), (f,))
*(c, f::CompositeFunction) = CompositeFunction(map(x -> c .* x, f.cs), f.fs)
+(f₁::GeoStatsFunction, f₂::GeoStatsFunction) = CompositeFunction((1, 1), (f₁, f₂))
+(f₁::CompositeFunction, f₂::GeoStatsFunction) = CompositeFunction((f₁.cs..., 1), (f₁.fs..., f₂))
+(f₁::GeoStatsFunction, f₂::CompositeFunction) = CompositeFunction((1, f₂.cs...), (f₁, f₂.fs...))
+(f₁::CompositeFunction, f₂::CompositeFunction) = CompositeFunction((f₁.cs..., f₂.cs...), (f₁.fs..., f₂.fs...))

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, f::CompositeFunction)
  O = IOContext(io, :compact => true)
  coeffs = 1 .* raw.(f.cs)
  models = nameof.(typeof.(f.fs))
  lines = ["$c × $f" for (c, f) in zip(coeffs, models)]
  print(O, join(lines, " + "))
end

function Base.show(io::IO, ::MIME"text/plain", f::CompositeFunction)
  O = IOContext(io, :compact => true)
  coeffs = 1 .* raw.(f.cs)
  models = f.fs
  name = "CompositeFunction"
  header = isisotropic(f) ? name : name * " (anisotropic)"
  println(O, header)
  println(O, "  structures")
  lines = ["    └─$f" for f in models]
  print(O, join(lines, "\n"))
  println(O)
  println(O, "  coefficients")
  lines = ["    └─$c" for c in coeffs]
  print(O, join(lines, "\n"))
end
