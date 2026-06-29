# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalVariogramSurface(geotable, vars;
                              normal=Vec(0,0,1), nangs=50,
                              ptol=0.5u"m", dtol=0.5u"m",
                              [options])

Given a `normal` direction, estimate the (cross-)variogram of variables `vars`
stored in `geotable` along all directions in the corresponding plane of variation.

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
  vars;
  normal=Vec(0, 0, 1),
  nangs=50,
  ptol=0.5u"m",
  dtol=0.5u"m",
  kwargs...
)
  # sanity checks
  @assert nangs > 1 "nangs must be greater than one"

  # deterministic results
  rng = Xoshiro(123)

  Dim = embeddim(domain(data))

  # basis for surface
  if Dim == 2
    planes = [data]
    u, v = Vec(1.0, 0.0), Vec(0.0, 1.0)
  elseif Dim == 3
    subset = partition(rng, data, PlanePartition(normal; tol=ptol))
    planes = collect(subset)
    u, v = Meshes.householderbasis(normal)
  else
    throw(ArgumentError("variogram surface only supported in 2D or 3D"))
  end

  # polar angles for half plane (variogram is symmetric)
  θs = collect(range(0, stop=π, length=nangs))

  # estimate directional variograms across planes
  gs = map(θs) do θ
    dir = DirectionPartition(cos(θ) * u + sin(θ) * v; tol=dtol)
    g(plane) = EmpiricalVariogram(partition(rng, plane, dir), vars; kwargs...)
    tmapreduce(g, merge, planes)
  end

  # polar radii
  rs = first(gs).abscissas

  # surface values
  zs = [g.ordinates for g in gs]

  EmpiricalVariogramSurface(θs, rs, zs)
end

nvariates(g::EmpiricalVariogramSurface) = size(first(g.zs), 1)

issymmetric(::Type{<:EmpiricalVariogramSurface}) = true
