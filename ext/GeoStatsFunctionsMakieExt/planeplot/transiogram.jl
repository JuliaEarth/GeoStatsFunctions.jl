# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _planeplot(
  t::EmpiricalTransioplane;

  # common options
  colormap=:viridis,
  maxlag=nothing,

  # transiogram options
  levels=nothing
)
  # polar angle
  θs = t.θs

  # polar radius
  rs = t.rs

  # hide hole at center
  rs = [zero(eltype(rs)); rs]

  # transioplane values
  zs = t.zs

  # number of labels
  L = size(first(zs), 1)

  # retrieve labels
  l = isnothing(levels) ? (1:L) : levels

  fig = Makie.Figure()
  for i in 1:L, j in 1:L
    lᵢ, lⱼ = l[i], l[j]
    ax = Makie.PolarAxis(fig[i, j], title="$lᵢ → $lⱼ")

    # values in matrix form
    Z = reduce(hcat, getindex.(zs, i, j))

    # hide hole at center
    Z = [Z[1:1, :]; Z]

    # transpose for plotting
    Z = transpose(Z)

    Makie.surface!(ax, θs, rs, Z, colormap=colormap, shading=Makie.NoShading)
  end

  fig
end

function planeplot(
  t::Transiogram;

  # common options
  colormap=:viridis,
  maxlag=nothing,

  # geometric options
  normal=nothing,
  nlags=20,
  nangs=50,

  # transiogram options
  levels=nothing
)
  # auxiliary parameters
  b = metricball(t)
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
  H = isnothing(maxlag) ? _maxlag(t) : _addunit(maxlag, u"m")

  # polar angles
  θs = range(0, stop=2π, length=2nangs)

  # polar radius
  rs = range(1e-6 * oneunit(H), stop=H, length=nlags)

  # direction vector as a function of polar angle
  dir(θ) = cos(θ) * u⃗ + sin(θ) * v⃗

  Z = [t(p, p + ustrip(r) * dir(θ)) for θ in θs, r in rs]

  # hide hole at center
  rs = [zero(eltype(rs)); rs]
  Z = [Z[:, 1:1] Z]

  # number of labels
  L = size(first(Z), 1)

  # retrieve labels
  l = isnothing(levels) ? (1:L) : levels

  fig = Makie.Figure()
  for i in 1:L, j in 1:L
    lᵢ, lⱼ = l[i], l[j]
    ax = Makie.PolarAxis(fig[i, j], title="$lᵢ → $lⱼ")

    # values in matrix form
    Zᵢⱼ = getindex.(Z, i, j)

    Makie.surface!(ax, θs, rs, Zᵢⱼ, colormap=colormap, shading=Makie.NoShading)
  end

  fig
end
