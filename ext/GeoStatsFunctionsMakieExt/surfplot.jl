# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function surfplot(gp::GridPos, f; labels=nothing, kwargs...)
  # variable names
  n = nvariates(f)
  l = isnothing(labels) ? (1:n) : labels

  layout = get_layout(gp)

  # initialize polar axis grid
  for i in 1:n, j in 1:n
    title = n > 1 ? "$(l[i]) → $(l[j])" : ""
    Makie.PolarAxis(layout[i, j], title=title)
  end

  # fill figure with plots
  surfplot!(gp, f; kwargs...)

  gp
end

function surfplot(f; figure=(;), kwargs...)
  fig = Makie.Figure(; figure...)
  surfplot(fig, f; kwargs...)
  fig
end

function surfplot!(
  gp::GridPos,
  f::GeoStatsFunction;
  # common options
  colormap=:viridis,
  maxlag=nothing,
  # geometric options
  normal=nothing,
  nlags=20,
  nangs=50
)
  # auxiliary parameters
  n = nvariates(f)
  m = _ncoords(f)
  U = typeof(range(f))

  # basis for plane
  if isnothing(normal)
    if m == 3
      d = 3
      u⃗ = Vec(U(1), U(0), U(0))
      v⃗ = Vec(U(0), U(1), U(0))
    else
      d = 2
      u⃗ = Vec(U(1), U(0))
      v⃗ = Vec(U(0), U(1))
    end
  else
    n⃗ = Vec(normal)
    d = length(n⃗)
    @assert d == m "normal vector is not compatible with given function"
    @assert d == 3 "normal vector must have 3 coordinates"
    u⃗, v⃗ = Meshes.householderbasis(n⃗)
  end

  # reference point
  p = Point(ntuple(i -> U(0), d))

  # direction vector as a function of polar angle
  dir(θ) = cos(θ) * u⃗ + sin(θ) * v⃗

  # maximum lag
  hmax = isnothing(maxlag) ? _maxlag(f) : GeoStatsFunctions.aslen(maxlag)

  # polar radius
  rs = range(1e-6 * oneunit(hmax), stop=hmax, length=nlags)

  # polar angle
  θs = if issymmetric(f)
    range(0, stop=π, length=nangs)
  else
    range(0, stop=2π, length=2nangs)
  end

  # evaluate function
  Z = [f(p, p + ustrip(r) * dir(θ)) for θ in θs, r in rs]

  # exploit symmetry
  if issymmetric(f)
    θs = range(0, 2π, length=2length(θs))
    Z = [Z; Z]
  end

  # hide hole at center
  rs = [zero(eltype(rs)); rs]
  Z = [Z[:, 1:1] Z]

  layout = get_layout(gp)
  for i in 1:n, j in 1:n
    ax = Makie.content(layout[i, j])
    Zᵢⱼ = getindex.(Z, i, j)
    Makie.surface!(ax, θs, ustrip.(rs), ustrip.(Zᵢⱼ), colormap=colormap, shading=false)
  end
  gp
end

_ncoords(f) = length(radii(metricball(f)))
_ncoords(::CarleTransiogram{N}) where {N} = N

function surfplot!(
  gp::GridPos,
  f::EmpiricalGeoStatsSurface;
  # common options
  colormap=:viridis,
  maxlag=nothing
)
  # auxiliary parameters
  n = nvariates(f)

  # polar radius
  rs = f.rs

  # polar angle
  θs = f.θs

  # function values
  zs = f.zs

  # maximum lag
  hmax = isnothing(maxlag) ? _maxlag(f) : GeoStatsFunctions.aslen(maxlag)

  # exploit symmetry
  if issymmetric(f)
    θs = range(0, stop=2π, length=2length(θs))
    zs = [zs; zs]
  end

  # discard above maximum lag
  is = findall(≤(hmax), rs)
  rs = rs[is]

  # hide hole at center
  rs = [zero(eltype(rs)); rs]

  layout = get_layout(gp)
  for i in 1:n, j in 1:n
    ax = Makie.content(layout[i, j])

    # values in matrix form
    Zᵢⱼ = _istransiogram(f) ? getindex.(zs, i, j) : zs
    Z = reduce(hcat, Zᵢⱼ)

    # discard above maximum lag
    Z = Z[is, :]

    # hide hole at center
    Z = [Z[1:1, :]; Z]

    # transpose for plotting
    Z = transpose(Z)

    Makie.surface!(ax, θs, ustrip.(rs), ustrip.(Z), colormap=colormap, shading=false)
  end
  gp
end
