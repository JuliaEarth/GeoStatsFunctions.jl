# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalVariogramSurface(data, var₁, var₂=var₁;
                              normal=Vec(0,0,1), nangs=50,
                              ptol=0.5u"m", dtol=0.5u"m",
                              [options])

Given a `normal` direction, estimate the (cross-)variogram of variables
`var₁` and `var₂` along all directions in the corresponding plane of variation.

Optionally, specify the tolerance `ptol` in length units for the plane partition,
the tolerance `dtol` in length units for the direction partition, the number of
angles `nangs` in the plane, and forward the `options` to the underlying
[`EmpiricalVariogram`](@ref).
"""
struct EmpiricalVariogramSurface{T,R,Z} <: EmpiricalGeoStatsSurface
  θs::Vector{T}
  rs::Vector{R}
  zs::Vector{Z}
end

function EmpiricalVariogramSurface(
  data::AbstractGeoTable,
  var₁,
  var₂=var₁;
  normal=Vec(0, 0, 1),
  nangs=50,
  ptol=0.5u"m",
  dtol=0.5u"m",
  kwargs...
)
  # sanity checks
  @assert nangs > 1 "nangs must be greater than one"

  # deterministic results
  rng = MersenneTwister(123)

  Dim = embeddim(domain(data))

  # basis for surface
  if Dim == 2
    planes = [data]
    u, v = Vec(1.0, 0.0), Vec(0.0, 1.0)
  elseif Dim == 3
    subset = partition(rng, data, PlanePartition(normal, tol=ptol))
    planes = collect(subset)
    u, v = householderbasis(normal)
  else
    throw(ArgumentError("variogram surface only supported in 2D or 3D"))
  end

  # polar angles for half plane (variogram is symmetric)
  θs = collect(range(0, stop=π, length=nangs))

  # estimate directional variograms across planes
  γs = map(θs) do θ
    dir = DirectionPartition(cos(θ) * u + sin(θ) * v, tol=dtol)
    γ(plane) = EmpiricalVariogram(partition(rng, plane, dir), var₁, var₂; kwargs...)
    tmapreduce(γ, merge, planes)
  end

  # polar radii
  rs = first(γs).abscissas

  # surface values
  zs = [γ.ordinates for γ in γs]

  EmpiricalVariogramSurface(θs, rs, zs)
end

nvariates(::Type{<:EmpiricalVariogramSurface}) = 1

issymmetric(::Type{<:EmpiricalVariogramSurface}) = true
