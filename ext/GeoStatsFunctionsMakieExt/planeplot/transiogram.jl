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
