# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function planeplot(g::EmpiricalVarioplane; colormap=:viridis)
  # polar angle
  θs = g.θs

  # polar radius
  rs = g.rs

  # varioplane values
  hs = g.hs

  # values in matrix form
  H = reduce(hcat, hs)

  # exploit symmetry
  θs = range(0, 2π, length=2 * length(θs))
  H = [H H]

  # hide hole at center
  rs = [zero(eltype(rs)); rs]
  H = [H[1:1, :]; H]

  # transpose for plotting
  H = transpose(H)

  fig = Makie.Figure()
  ax = Makie.PolarAxis(fig[1, 1])
  Makie.surface!(ax, θs, rs, H, colormap=colormap, shading=Makie.NoShading)
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
  hs = t.hs

  # number of labels
  L = size(first(hs), 1)

  # retrieve labels
  l = isnothing(levels) ? (1:L) : levels

  fig = Makie.Figure()
  for i in 1:L, j in 1:L
    lᵢ, lⱼ = l[i], l[j]
    ax = Makie.PolarAxis(fig[i, j], title="$lᵢ → $lⱼ")

    # values in matrix form
    h = getindex.(hs, i, j)
    H = reduce(hcat, h)

    # hide hole at center
    H = [H[1:1, :]; H]

    # transpose for plotting
    H = transpose(H)

    Makie.surface!(ax, θs, rs, H, colormap=colormap, shading=Makie.NoShading)
  end

  fig
end
