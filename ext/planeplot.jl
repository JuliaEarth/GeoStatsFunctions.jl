# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function planeplot(g::EmpiricalVarioplane; colormap=:viridis)
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

function planeplot(t::EmpiricalTransioplane; colormap=:viridis, levels=nothing)
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

function planeplot(γ::Variogram;
  # geometric options
  normal=Vec(0, 0, 1),
  nlags=20, nangs=50,
  maxlag=nothing,
  # austhetic options
  colormap=:viridis
)
  # basis for varioplane
  u, v = householderbasis(normal)

  # maximum lag
  H = if isnothing(maxlag)
    _maxlag(γ)
  else
    _addunit(maxlag, u"m")
  end

  # polar angles
  θs = range(0, stop=π, length=nangs)

  # polar radius
  rs = range(1e-6, stop=ustrip(H), length=nlags)

  # direction vector as a function of polar angle
  dir(θ) = cos(θ) * u + sin(θ) * v

  O = Point(0, 0, 0)
  G = [γ(O, O + r * dir(θ)) for θ in θs, r in rs]

  # exploit symmetry
  θs = range(0, 2π, length=2 * length(θs))
  G = [G; G]

  # hide hole at center
  rs = [zero(eltype(rs)); rs]
  G = [G[:, 1:1] G]

  fig = Makie.Figure()
  ax = Makie.PolarAxis(fig[1, 1])
  Makie.surface!(ax, θs, rs, G, colormap=colormap, shading=Makie.NoShading)
  fig
end
