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
    isisotropic(f)

Tells whether or not the geostatistical function `f` is isotropic.
"""
isisotropic(f::GeoStatsFunction) = isisotropic(metricball(f))

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
    returntype(f, g₁, g₂)

Return type of f(g₁, g₂).
"""
returntype(f::GeoStatsFunction, g₁, g₂) = typeof(f(g₁, g₂))

"""
    pairwise(f, domain)

Evaluate geostatistical function `f` between all elements in the `domain`.
"""
function pairwise(f::GeoStatsFunction, domain)
  g = first(domain)
  n = length(domain)
  T = returntype(f, g, g)
  F = Matrix{T}(undef, n, n)
  pairwise!(F, f, domain)
end

"""
    pairwise!(F, f, domain)

Evaluates geostatistical function `f` between all elements in the `domain` in-place, filling the matrix `F`.
"""
function pairwise!(F, f::GeoStatsFunction, domain)
  n = length(domain)
  @inbounds for j in 1:n
    gⱼ = domain[j]
    sⱼ = _sample(f, gⱼ)
    for i in (j + 1):n
      gᵢ = domain[i]
      sᵢ = _sample(f, gᵢ)
      F[i, j] = mean(f(pᵢ, pⱼ) for pᵢ in sᵢ, pⱼ in sⱼ)
    end
    F[j, j] = mean(f(pⱼ, pⱼ) for pⱼ in sⱼ, pⱼ in sⱼ)
    for i in 1:(j - 1)
      F[i, j] = F[j, i] # leverage the symmetry
    end
  end
  F
end

"""
    pairwise(f, domain₁, domain₂)

Evaluate geostatistical function `f` between all elements of `domain₁` and `domain₂`.
"""
function pairwise(f::GeoStatsFunction, domain₁, domain₂)
  g₁ = first(domain₁)
  g₂ = first(domain₂)
  m = length(domain₁)
  n = length(domain₂)
  T = returntype(f, g₁, g₂)
  F = Matrix{T}(undef, m, n)
  pairwise!(F, f, domain₁, domain₂)
end

"""
    pairwise!(F, f, domain₁, domain₂)

Evaluates geostatistical function `f` between all elements of `domain₁` and `domain₂` in-place, filling the matrix `F`.
"""
function pairwise!(F, f::GeoStatsFunction, domain₁, domain₂)
  m = length(domain₁)
  n = length(domain₂)
  @inbounds for j in 1:n
    gⱼ = domain₂[j]
    sⱼ = _sample(f, gⱼ)
    for i in 1:m
      gᵢ = domain₁[i]
      sᵢ = _sample(f, gᵢ)
      F[i, j] = mean(f(pᵢ, pⱼ) for pᵢ in sᵢ, pⱼ in sⱼ)
    end
  end
  F
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("theoretical/variogram.jl")
include("theoretical/covariance.jl")
include("theoretical/transiogram.jl")

# -----------------
# HELPER FUNCTIONS
# -----------------

_sample(::GeoStatsFunction, p::Point) = [p]

function _sample(f::GeoStatsFunction, g::Geometry)
  α = _spacing(f, g)
  rng = MersenneTwister(123)
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
