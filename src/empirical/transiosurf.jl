# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalTransiogramSurface(θs, rs, zs, vs)

Type that stores all information about the empirical (a.k.a. experimental)
transiogram surface. It is returned by the [`transiogramsurface`](@ref) function
and is considered an internal type (i.e., not intended for end-users).
"""
struct EmpiricalTransiogramSurface{T,R,Z} <: EmpiricalGeoStatsSurface
  θs::Vector{T}
  rs::Vector{R}
  zs::Vector{Z}
  vs::Vector{Symbol}
end

issymmetric(::Type{<:EmpiricalTransiogramSurface}) = false

nvariables(t::EmpiricalTransiogramSurface) = size(first(t.zs), 1)

variables(t::EmpiricalTransiogramSurface) = t.vs

# -------------------
# END-USER INTERFACE
# -------------------

"""
    transiogramsurface(geotable, [var];
                       normal=Vec(0,0,1), nangs=50,
                       planetol=0.5u"m", dirtol=0.5u"m",
                       [options])

Given a `normal` direction, estimate the transiogram of categorical variable `var`
stored in `geotable` along all directions in the corresponding plane of variation.

Optionally, specify the tolerance `planetol` in length units for the plane partition,
the tolerance `dirtol` in length units for the direction partition, the number of
angles `nangs` in the plane, and forward the `options` to the underlying
[`transiogram`](@ref).
"""
function transiogramsurface(
  data::AbstractGeoTable,
  var=1;
  normal=Vec(0, 0, 1),
  nangs=50,
  planetol=0.5u"m",
  dirtol=0.5u"m",
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
    subset = partition(rng, data, PlanePartition(normal; tol=planetol))
    planes = collect(subset)
    u, v = Meshes.householderbasis(normal)
  else
    throw(ArgumentError("transiogram surface only supported in 2D or 3D"))
  end

  # polar angles for full plane (transiogram is asymmetric)
  θs = collect(range(0, stop=2π, length=2nangs))

  # estimate directional transiograms across planes
  ts = map(θs) do θ
    dir = DirectionPartition(cos(θ) * u + sin(θ) * v; tol=dirtol)
    t(plane) = transiogram(partition(rng, plane, dir), var; kwargs...)
    tmapreduce(t, merge, planes)
  end

  # polar radii
  rs = first(ts).abscissas

  # surface values
  zs = [t.ordinates for t in ts]

  # variable names
  vs = variables(first(ts))

  EmpiricalTransiogramSurface(θs, rs, zs, vs)
end
