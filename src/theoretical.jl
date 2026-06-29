# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeoStatsFunction

A theoretical geostatistical function (e.g. variogram, covariance).
"""
abstract type GeoStatsFunction end

"""
    constructor(f)

Return the type constructor of the geostatistical function `f`.
"""
function constructor end

"""
    isstationary(f)

Check if geostatistical function `f` possesses the 2nd-order stationary property.
"""
isstationary(f::GeoStatsFunction) = isstationary(typeof(f))

"""
    isisotropic(f)

Tells whether or not the geostatistical function `f` is isotropic.
"""
isisotropic(f::GeoStatsFunction) = hasequalradii(metricball(f))

"""
    issymmetric(f)

Tell whether or not the geostatistical function `f` is symmetric.
"""
issymmetric(f::GeoStatsFunction) = issymmetric(typeof(f))

"""
    isbanded(f)

Tells whether or not the geostatistical function `f` produces a banded matrix.
"""
isbanded(f::GeoStatsFunction) = isbanded(typeof(f))

"""
    metricball(f)

Return the metric ball of the geostatistical function `f`.
"""
metricball(f::GeoStatsFunction) = f.ball

"""
    range(f)

Return the maximum effective range of the geostatistical function `f`.
"""
Base.range(f::GeoStatsFunction) = maximum(radii(metricball(f)))

"""
    scale(f, s)

Scale the metric ball of the geostatistical function `f`
with a strictly positive scaling factor `s`.
"""
function scale end

"""
    nvariables(f)

Return the number of (co)variables of the geostatistical function `f`.
"""
nvariables(f::GeoStatsFunction) = nvariables(typeof(f))

"""
    variables(f)

Return the names of the (co)variables of the geostatistical function `f`.
"""
variables(f::GeoStatsFunction) = defaultvariables(nvariables(f))

"""
    f(p₁, p₂)

Evaluate the geostatistical function at points `p₁` and `p₂`.
"""
function (f::GeoStatsFunction)(p₁::Point, p₂::Point)
  d = metric(metricball(f))
  h = evaluate(d, p₁, p₂)
  f(h)
end

"""
    f(g₁, p₂)

Evaluate the geostatistical function at geometry `g₁` and point `p₂`.
"""
function (f::GeoStatsFunction)(g₁::Geometry, p₂::Point)
  s₁ = _sample(f, g₁)
  mean(f(p₁, p₂) for p₁ in s₁)
end

"""
    f(p₁, g₂)

Evaluate the geostatistical function at point `p₁` and geometry `g₂`.
"""
(f::GeoStatsFunction)(p₁::Point, g₂::Geometry) = f(g₂, p₁)

"""
    f(g₁, g₂)

Evaluate the geostatistical function at geometries `g₁` and `g₂`.
"""
function (f::GeoStatsFunction)(g₁::Geometry, g₂::Geometry)
  s₁ = _sample(f, g₁)
  s₂ = _sample(f, g₂)
  mean(f(p₁, p₂) for p₁ in s₁, p₂ in s₂)
end

"""
    pairwise(f, domain)

Evaluate geostatistical function `f` between all elements in the `domain`,
filling a matrix with unitless values for use in linear systems.
"""
function pairwise(f::GeoStatsFunction, domain)
  T, (m, n, k) = matrixparams(f, domain)
  F = Matrix{T}(undef, (m * k, n * k))
  pairwise!(F, f, domain)
end

"""
    pairwise(f, domain₁, domain₂)

Evaluate geostatistical function `f` between all elements in `domain₁` and `domain₂`
filling a matrix with unitless values for use in linear systems.
"""
function pairwise(f::GeoStatsFunction, domain₁, domain₂)
  T, (m, n, k) = matrixparams(f, domain₁, domain₂)
  F = Matrix{T}(undef, (m * k, n * k))
  pairwise!(F, f, domain₁, domain₂)
end

"""
    pairwise!(F, f, domain)

Evaluates geostatistical function `f` between all elements in the `domain` in-place,
filling the matrix `F` with unitless values for use in linear systems.
"""
pairwise!(F, f::GeoStatsFunction, domain) = pairwise!(F, f, domain, domain)

"""
    pairwise!(F, f, domain₁, domain₂)

Evaluates geostatistical function `f` between all elements of `domain₁` and `domain₂` in-place,
filling the matrix `F` with unitless values for use in linear systems.
"""
function pairwise!(F, f::GeoStatsFunction, domain₁, domain₂)
  _, (m, n, k) = matrixparams(f, domain₁, domain₂)
  @inbounds for j in 1:n
    gⱼ = domain₂[j]
    sⱼ = _sample(f, gⱼ)
    for i in 1:m
      gᵢ = domain₁[i]
      sᵢ = _sample(f, gᵢ)
      Fᵢⱼ = ustrip.(mean(f(pᵢ, pⱼ) for pᵢ in sᵢ, pⱼ in sⱼ))
      F[((i - 1) * k + 1):(i * k), ((j - 1) * k + 1):(j * k)] .= Fᵢⱼ
    end
  end
  F
end

"""
    matrixparams(f, domain)

Return the parameters used to assemble the matrix for the
geostatistical function `f` between all elements in the `domain`.
"""
matrixparams(f::GeoStatsFunction, domain) = matrixparams(f, domain, domain)

"""
    matrixparams(f, domain₁, domain₂)

Return the parameters used to assemble the matrix for the
geostatistical function `f` between all elements in `domain₁`
and `domain₂`.
"""
function matrixparams(f::GeoStatsFunction, domain₁, domain₂)
  m = length(domain₁)
  n = length(domain₂)
  k = nvariables(f)
  g₁ = first(domain₁)
  g₂ = first(domain₂)
  F = ustrip.(f(g₁, g₂))
  V = eltype(F)
  V, (m, n, k)
end

# -----------
# IO METHODS
# -----------

Base.summary(io::IO, f::GeoStatsFunction) = print(io, nameof(typeof(f)))

function Base.show(io::IO, f::GeoStatsFunction)
  ioctx = IOContext(io, :compact => true)
  summary(ioctx, f)
  print(ioctx, "(")
  _printfields(ioctx, f, singleline=true)
  print(ioctx, ")")
end

function Base.show(io::IO, ::MIME"text/plain", f::GeoStatsFunction)
  ioctx = IOContext(io, :compact => true, :limit => true)
  summary(ioctx, f)
  println(ioctx)
  _printfields(ioctx, f)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("theoretical/variogram.jl")
include("theoretical/covariance.jl")
include("theoretical/transiogram.jl")
include("theoretical/composite.jl")

# ---------------------
# THEORETICAL MATRICES
# ---------------------

include("theoretical/matrices.jl")

# -----------------
# HELPER FUNCTIONS
# -----------------

function _printfields(io, f; singleline=false)
  fields = fieldnames(typeof(f))
  len = length(fields)
  precommon, prelast = singleline ? ("", "") : ("├─ ", "└─ ")
  poscommon, poslast = singleline ? (", ", "") : ("\n", "")
  for (i, field) in enumerate(fields)
    pre = i == len ? prelast : precommon
    pos = i == len ? poslast : poscommon
    value = getfield(f, i)
    if value isa MetricBall
      print(io, pre)
      if hasequalradii(value)
        print(io, "range: ", first(radii(value)))
      else
        print(io, "ranges: ", radii(value))
        print(io, pos)
        print(io, pre)
        print(io, "rotation: ", rotation(value))
      end
      print(io, pos)
    else
      print(io, pre)
      print(io, field, ": ", value)
      print(io, pos)
    end
  end
end

_sample(::GeoStatsFunction, p::Point) = (p,)

function _sample(f::GeoStatsFunction, g::Geometry)
  rng = MersenneTwister(123)
  α = _spacing(f, g)
  sample(rng, g, MinDistanceSampling(α))
end

function _sample(f::GeoStatsFunction, s::Segment)
  s′ = Segment(s(0.05), s(0.95))
  n = _dims(f, s)
  sample(s′, RegularSampling(n...))
end

function _sample(f::GeoStatsFunction, q::Quadrangle)
  q′ = Quadrangle(q(0.05, 0.05), q(0.95, 0.05), q(0.95, 0.95), q(0.05, 0.95))
  n = _dims(f, q)
  sample(q′, RegularSampling(n...))
end

function _sample(f::GeoStatsFunction, h::Hexahedron)
  h′ = Hexahedron(
    h(0.05, 0.05, 0.05),
    h(0.95, 0.05, 0.05),
    h(0.95, 0.95, 0.05),
    h(0.05, 0.95, 0.05),
    h(0.05, 0.05, 0.95),
    h(0.95, 0.05, 0.95),
    h(0.95, 0.95, 0.95),
    h(0.05, 0.95, 0.95)
  )
  n = _dims(f, h)
  sample(h′, RegularSampling(n...))
end

function _spacing(f::GeoStatsFunction, g::Geometry)
  r = range(f)
  s = sides(boundingbox(g))
  l = minimum(filter(sᵢ -> sᵢ > zero(sᵢ), s))
  r > zero(r) ? min(r, l) / 3 : l / 3
end

function _dims(f::GeoStatsFunction, g::Geometry)
  s = sides(boundingbox(g))
  α = _spacing(f, g)
  n = ceil.(Int, s ./ α)
  ntuple(i -> iszero(n[i]) ? 1 : n[i], length(n))
end
