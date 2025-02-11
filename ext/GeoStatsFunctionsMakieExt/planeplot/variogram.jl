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
