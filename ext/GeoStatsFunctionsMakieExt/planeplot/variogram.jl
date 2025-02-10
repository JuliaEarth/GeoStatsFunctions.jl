# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _planeplot(
  g::EmpiricalVarioplane;

  # common options
  colormap=:viridis,
  maxlag=nothing
)
  # polar angle
  θs = g.θs

  # polar radius
  rs = g.rs

  # varioplane values
  zs = g.zs

  # values in matrix form
  Z = reduce(hcat, zs)

  # exploit symmetry
  θs = range(0, 2π, length=2 * length(θs))
  Z = [Z Z]

  # hide hole at center
  rs = [zero(eltype(rs)); rs]
  Z = [Z[1:1, :]; Z]

  # transpose for plotting
  Z = transpose(Z)

  fig = Makie.Figure()
  ax = Makie.PolarAxis(fig[1, 1])
  Makie.surface!(ax, θs, rs, Z, colormap=colormap, shading=Makie.NoShading)

  fig
end

function _planeplot(
  γ::Variogram;

  # common options
  colormap=:viridis,
  maxlag=nothing,

  # geometric options
  normal=nothing,
  nlags=20,
  nangs=50
)
  # auxiliary parameters
  b = metricball(γ)
  r = radii(b)
  n = length(r)
  U = eltype(r)

  # retrieve normal direction
  n⃗ = isnothing(normal) ? Vec(ntuple(i -> U(i == n), n)) : normal

  # sanity checks
  @assert length(n⃗) == n "normal direction must have $n coordinates"

  # basis for varioplane
  u⃗, v⃗ = if n == 3
    householderbasis(n⃗)
  else
    Vec(U(1), U(0)), Vec(U(0), U(1))
  end

  # reference point
  p = Point(ntuple(i -> U(0), n))

  # retrieve maximum lag
  H = isnothing(maxlag) ? _maxlag(γ) : _addunit(maxlag, u"m")

  # polar angles
  θs = range(0, stop=π, length=nangs)

  # polar radius
  rs = range(1e-6 * oneunit(H), stop=H, length=nlags)

  # direction vector as a function of polar angle
  dir(θ) = cos(θ) * u⃗ + sin(θ) * v⃗

  Z = [γ(p, p + ustrip(r) * dir(θ)) for θ in θs, r in rs]

  # exploit symmetry
  θs = range(0, 2π, length=2 * length(θs))
  Z = [Z; Z]

  # hide hole at center
  rs = [zero(eltype(rs)); rs]
  Z = [Z[:, 1:1] Z]

  fig = Makie.Figure()
  ax = Makie.PolarAxis(fig[1, 1])
  Makie.surface!(ax, θs, rs, Z, colormap=colormap, shading=Makie.NoShading)

  fig
end
