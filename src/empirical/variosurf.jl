# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalVariogramSurface(θs, rs, zs, vs)

Type that stores all information about the empirical (a.k.a. experimental)
variogram surface. It is returned by the [`variogramsurface`](@ref) function
and is considered an internal type (i.e., not intended for end-users).
"""
struct EmpiricalVariogramSurface{T,R,Z} <: EmpiricalGeoStatsSurface
  θs::Vector{T}
  rs::Vector{R}
  zs::Vector{Z}
  vs::Vector{Symbol}
end

issymmetric(::Type{<:EmpiricalVariogramSurface}) = true

nvariables(g::EmpiricalVariogramSurface) = size(first(g.zs), 1)

variables(g::EmpiricalVariogramSurface) = g.vs

# -------------------
# END-USER INTERFACE
# -------------------

"""
    variogramsurface(geotable, [vars];
                     normal=Vec(0,0,1), nangs=50,
                     ptol=0.5u"m", dtol=0.5u"m",
                     [options])

Given a `normal` direction, estimate the (cross-)variogram of variables `vars`
stored in `geotable` along all directions in the corresponding plane of variation.

Optionally, specify the tolerance `ptol` in length units for the plane partition,
the tolerance `dtol` in length units for the direction partition, the number of
angles `nangs` in the plane, and forward the `options` to the underlying
[`variogram`](@ref) calls.
"""
function variogramsurface(
  data::AbstractGeoTable,
  vars=1:(ncol(data) - 1);
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

  dim = embeddim(domain(data))

  # basis for surface
  if dim == 2
    planes = [data]
    u, v = Vec(1.0, 0.0), Vec(0.0, 1.0)
  elseif dim == 3
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
    g(plane) = variogram(partition(rng, plane, dir), vars; kwargs...)
    tmapreduce(g, merge, planes)
  end

  # polar radii
  rs = first(gs).abscissas

  # surface values
  zs = [g.ordinates for g in gs]

  # variable names
  vs = variables(first(gs))

  EmpiricalVariogramSurface(θs, rs, zs, vs)
end
