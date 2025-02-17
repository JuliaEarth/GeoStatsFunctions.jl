# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CompositeFunction(cs, fs)

A composite geostatistical function `f = c₁f₁ + c₂f₂ + ⋯ + cₙfₙ` with
coefficients `cs = (c₁, c₂, ..., cₙ)` and geostatistical functions
`fs = (f₁, f₂, ..., fₙ)`.
"""
struct CompositeFunction{CS,FS} <: GeoStatsFunction
  cs::CS
  fs::FS

  function CompositeFunction(cs::CS, fs::FS) where {CS,FS}
    # convert arrays to static arrays of floats
    cs′ = map(c -> _static(float(c)), cs)
    CS′ = typeof(cs′)
    new{CS′,FS}(cs′, fs)
  end
end

isstationary(f::CompositeFunction) = all(isstationary, f.fs)

isisotropic(f::CompositeFunction) = all(isisotropic, f.fs)

issymmetric(f::CompositeFunction) = all(issymmetric, f.cs) && all(issymmetric, f.fs)

isbanded(f::CompositeFunction) = all(isbanded, f.fs)

sill(f::CompositeFunction) = sum(f.cs .* map(sill, f.fs))

nugget(f::CompositeFunction) = sum(f.cs .* map(nugget, f.fs))

metricball(f::CompositeFunction) = metricball(argmax(fᵢ -> range(fᵢ), f.fs))

Base.range(f::CompositeFunction) = maximum(range(fᵢ) for fᵢ in f.fs)

nvariates(f::CompositeFunction) = maximum(size(cᵢ, 1) for cᵢ in f.cs)

scale(f::CompositeFunction, s::Real) = CompositeFunction(f.cs, map(fᵢ -> scale(fᵢ, s), f.fs))

function structures(f::CompositeFunction)
  ks, fs = f.cs, f.fs

  # total nugget and contributions
  cₒ = sum(@. ks * nugget(fs))
  cs = @. ks * (sill(fs) - nugget(fs))

  # discard nugget effect terms
  inds = findall(fᵢ -> !(fᵢ isa NuggetEffect), fs)
  cs, fs = cs[inds], fs[inds]

  # adjust sill and nugget
  fs = map(f -> first(last(structures(f))), fs)

  # strip units from coefficients
  ucₒ = ustrip.(cₒ)
  ucs = map(c -> ustrip.(c), cs)

  ucₒ, ucs, fs
end

(f::CompositeFunction)(h) = sum(f.cs .* map(fᵢ -> fᵢ(h), f.fs))
(f::CompositeFunction)(u::Point, v::Point) = sum(f.cs .* map(fᵢ -> fᵢ(u, v), f.fs))

# algebraic structure
*(c, f::GeoStatsFunction) = CompositeFunction((c,), (f,))
*(c, f::CompositeFunction) = CompositeFunction(map(cᵢ -> c .* cᵢ, f.cs), f.fs)
+(f₁::GeoStatsFunction, f₂::GeoStatsFunction) = CompositeFunction((1, 1), (f₁, f₂))
+(f₁::CompositeFunction, f₂::GeoStatsFunction) = CompositeFunction((f₁.cs..., 1), (f₁.fs..., f₂))
+(f₁::GeoStatsFunction, f₂::CompositeFunction) = CompositeFunction((1, f₂.cs...), (f₁, f₂.fs...))
+(f₁::CompositeFunction, f₂::CompositeFunction) = CompositeFunction((f₁.cs..., f₂.cs...), (f₁.fs..., f₂.fs...))

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, f::CompositeFunction)
  O = IOContext(io, :compact => true)
  coeffs = f.cs
  models = nameof.(typeof.(f.fs))
  lines = ["$c × $f" for (c, f) in zip(coeffs, models)]
  print(O, join(lines, " + "))
end

function Base.show(io::IO, ::MIME"text/plain", f::CompositeFunction)
  O = IOContext(io, :compact => true)
  coeffs = f.cs
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

# -----------------
# HELPER FUNCTIONS
# -----------------

_static(c) = c
_static(c::AbstractMatrix) = SMatrix{size(c)...}(c)
