# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _funplot(
  γ::EmpiricalVariogram;
  # common options
  color=:teal,
  size=1.5,
  maxlag=nothing,
  # empirical options
  pointsize=12,
  showtext=true,
  textsize=12,
  showhist=true,
  histcolor=:teal
)
  # retrieve coordinates and counts
  x = γ.abscissas
  y = γ.ordinates
  n = γ.counts

  # retrieve maximum lag
  H = isnothing(maxlag) ? last(x) : _addunit(maxlag, u"m")

  # discard empty bins
  x = x[n .> 0]
  y = y[n .> 0]
  n = n[n .> 0]

  # discard above maximum lag
  n = n[x .≤ H]
  y = y[x .≤ H]
  x = x[x .≤ H]

  # initialize figure and axis
  fig = Makie.Figure()
  ax = Makie.Axis(fig[1, 1])

  # visualize frequencies as bars
  if showhist
    f = n * (maximum(y) / maximum(n)) / 10
    Makie.barplot!(ax, x, f, color=histcolor, alpha=0.3, gap=0.0)
  end

  # visualize variogram
  Makie.scatterlines!(ax, x, y, color=color, markersize=pointsize, linewidth=size)

  # visualize text counts
  if showtext
    text = string.(n)
    Makie.text!(ax, x, y, text=text, fontsize=textsize)
  end

  fig
end
