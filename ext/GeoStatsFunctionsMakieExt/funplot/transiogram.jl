# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _funplot(
  t::EmpiricalTransiogram;
  # common options
  color=:teal,
  size=1.5,
  maxlag=nothing,
  labels=nothing,
  # empirical options
  pointsize=12,
  showtext=true,
  textsize=12,
  showhist=true,
  histcolor=:teal
)
  # number of labels
  L = Base.size(t.ordinates, 1)

  # retrieve labels
  l = isnothing(labels) ? (1:L) : labels

  fig = Makie.Figure()
  for i in 1:L, j in 1:L
    lᵢ, lⱼ = l[i], l[j]
    ax = Makie.Axis(fig[i, j])

    # retrieve coordinates and counts
    x = t.abscissas
    y = t.ordinates[i, j]
    n = t.counts

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

    # visualize frequencies as bars
    if showhist
      f = n * (maximum(y) / maximum(n)) / 10
      Makie.barplot!(ax, x, f, color=histcolor, alpha=0.3, gap=0.0)
    end

    # visualize transiogram
    Makie.scatterlines!(ax, x, y, color=color, markersize=pointsize, linewidth=size, label="$lᵢ → $lⱼ")

    # visualize text counts
    if showtext
      text = string.(n)
      Makie.text!(ax, x, y, text=text, fontsize=textsize)
    end

    Makie.axislegend(position=i == j ? :rt : :rb)
  end
  fig
end
