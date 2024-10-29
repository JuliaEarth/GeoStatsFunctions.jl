# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalTransioplane(data, var;
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
struct EmpiricalTransioplane{T,V}
  θs::Vector{T}
  ts::Vector{V}
end

function EmpiricalTransioplane(
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

  # basis for transioplane
  if Dim == 2
    planes = [data]
    u, v = Vec(1.0, 0.0), Vec(0.0, 1.0)
  elseif Dim == 3
    subset = partition(rng, data, PlanePartition(normal, tol=ptol))
    planes = collect(subset)
    u, v = householderbasis(normal)
  else
    @error "transioplane only supported in 2D or 3D"
  end

  # loop over all angles of the plane
  θs = collect(range(0, stop=2π, length=nangs))
  ts = map(θs) do θ
    dir = DirectionPartition(cos(θ) * u + sin(θ) * v, tol=dtol)
    γ(plane) = EmpiricalTransiogram(partition(rng, plane, dir), var₁, var₂; kwargs...)
    tmapreduce(γ, merge, planes)
  end

  EmpiricalTransioplane(θs, ts)
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, ::EmpiricalTransioplane)
  print(io, "EmpiricalTransioplane")
end

function Base.show(io::IO, ::MIME"text/plain", γ::EmpiricalTransioplane)
  θs = [@sprintf "%.2f" rad2deg(θ) for θ in γ.θs]
  nθ = length(θs)
  lines = ["  └─$(θ)°" for θ in θs]
  lines = length(lines) > 11 ? vcat(lines[1:5], ["  ⋮"], lines[(end - 4):end]) : lines
  println(io, γ)
  println(io, "  $nθ angles")
  print(io, join(lines, "\n"))
end
