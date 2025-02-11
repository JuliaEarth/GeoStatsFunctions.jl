# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalTransiogramSurface(data, var;
                                normal=Vec(0,0,1), nangs=50,
                                ptol=0.5u"m", dtol=0.5u"m",
                                [parameters])

Given a `normal` direction, estimate the transiogram of variable `var`
along all directions in the corresponding plane of variation.

Optionally, specify the tolerance `ptol` in length units for the plane partition,
the tolerance `dtol` in length units for the direction partition, the number of
angles `nangs` in the plane, and forward the `parameters` to the underlying
[`EmpiricalTransiogram`](@ref).
"""
struct EmpiricalTransiogramSurface{T,R,Z} <: EmpiricalGeoStatsSurface
  θs::Vector{T}
  rs::Vector{R}
  zs::Vector{Z}
end

function EmpiricalTransiogramSurface(
  data::AbstractGeoTable,
  var;
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
    throw(ArgumentError("transiogram surface only supported in 2D or 3D"))
  end

  # polar angles for full plane (transiogram is asymmetric)
  θs = collect(range(0, stop=2π, length=2nangs))

  # estimate directional transiograms across planes
  ts = map(θs) do θ
    dir = DirectionPartition(cos(θ) * u + sin(θ) * v, tol=dtol)
    t(plane) = EmpiricalTransiogram(partition(rng, plane, dir), var; kwargs...)
    tmapreduce(t, merge, planes)
  end

  # polar radii
  rs = first(ts).abscissas

  # surface values
  zs = [t.ordinates for t in ts]

  EmpiricalTransiogramSurface(θs, rs, zs)
end

nvariates(t::EmpiricalTransiogramSurface) = size(first(t.zs), 1)

issymmetric(::Type{<:EmpiricalTransiogramSurface}) = false
