# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function surfplot(f; names=nothing, kwargs...)
  # initialize figure
  n = nvariates(f)
  v = isnothing(names) ? (1:n) : names
  fig = Makie.Figure()
  for i in 1:n, j in 1:n
    ax = Makie.PolarAxis(fig[i, j])
    ax.title = n > 1 ? "$(v[i]) → $(v[j])" : ""
  end

  # fill figure with plots
  surfplot!(fig, f; kwargs...)
end

function surfplot!(
  fig::Makie.Figure,
  f::GeoStatsFunction;
  # common options
  colormap=:viridis,
  maxlag=nothing,
  # geometric options
  normal=nothing,
  nlags=20,
  nangs=50
)
  # maximum lag
  hmax = isnothing(maxlag) ? _maxlag(f) : GeoStatsFunctions.aslen(maxlag)

  # basis for plane
  m = _ncoords(f)
  U = typeof(range(f))
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

  # add plots to axes
  n = nvariates(f)
  I = LinearIndices((n, n))
  for i in 1:n, j in 1:n
    ax = fig.content[I[j, i]]

    # values in matrix form
    Zᵢⱼ = getindex.(Z, i, j)

    # compute contour levels
    levels = quantile(Zᵢⱼ, range(0, 1, length=20))

    Makie.contour!(ax, θs, ustrip.(u"m", rs), ustrip.(Zᵢⱼ); colormap, levels, labels=true)
  end

  fig
end

_ncoords(f) = length(radii(metricball(f)))
_ncoords(::CarleTransiogram{N}) where {N} = N

function surfplot!(
  fig::Makie.Figure,
  f::EmpiricalGeoStatsSurface;
  # common options
  colormap=:viridis,
  maxlag=nothing
)
  # maximum lag
  hmax = isnothing(maxlag) ? _maxlag(f) : GeoStatsFunctions.aslen(maxlag)

  # polar radius
  rs = f.rs

  # polar angle
  θs = f.θs

  # function values
  zs = f.zs

  # exploit symmetry
  if issymmetric(f)
    θs = range(0, stop=2π, length=2length(θs))
    zs = [zs; zs]
  end

  # discard above maximum lag
  is = findall(≤(hmax), rs)
  rs = rs[is]

  # add plots to axes
  n = nvariates(f)
  I = LinearIndices((n, n))
  for i in 1:n, j in 1:n
    ax = fig.content[I[j, i]]

    # values in matrix form
    Zᵢⱼ = getindex.(zs, i, j)
    Z = reduce(hcat, Zᵢⱼ)

    # discard above maximum lag
    Z = Z[is, :]

    # transpose for plotting
    Z = transpose(Z)

    # compute contour levels
    levels = quantile(Z, range(0, 1, length=20))

    Makie.contour!(ax, θs, ustrip.(u"m", rs), ustrip.(Z); colormap, levels, labels=true)
  end

  fig
end
