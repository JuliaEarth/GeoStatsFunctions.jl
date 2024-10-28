# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ----------
# EMPIRICAL
# ----------

function transioplot(
  t::EmpiricalTransiogram;
  # common transiogram options
  color=:slategray,
  size=1.5,
  maxlag=nothing,
  levels=nothing,

  # empirical transiogram options
  pointsize=12,
  showtext=true,
  textsize=12,
  showhist=true,
  histcolor=:slategray
)
  # number of labels
  L = Base.size(t.ordinates, 1)

  # retrieve labels
  l = isnothing(levels) ? (1:L) : levels

  fig = Makie.Figure()
  for i in 1:L, j in 1:L
    lᵢ, lⱼ = l[i], l[j]
    ax = Makie.Axis(fig[i, j])

    # retrieve coordinates and counts
    x = t.abscissas
    y = t.ordinates[i, j]
    n = t.counts

    # discard empty bins
    x = x[n .> 0]
    y = y[n .> 0]
    n = n[n .> 0]

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

# ------------
# THEORETICAL
# ------------

function transioplot(
  t::Transiogram;
  # common transiogram options
  color=:slategray,
  size=1.5,
  maxlag=nothing,
  levels=nothing
)
  # retrieve maximum lag
  H = if isnothing(maxlag)
    _maxlag(t)
  else
    _addunit(maxlag, u"m")
  end

  # transiogram up to maximum lag
  hs = range(zero(H), stop=H, length=100)
  ts = t.(hs)

  # categorical labels
  L = Base.size(first(ts), 1)
  l = isnothing(levels) ? (1:L) : levels

  # effective ranges
  r = GeoStatsFunctions.ranges(t)

  # relative proportions at h → ∞
  p = normalize(diag(t(100H)), 1)

  # base transiogram model
  t̂ = MatrixExponentialTransiogram(r, p)
  t̂s = t̂.(hs)

  fig = Makie.Figure()
  for i in 1:L, j in 1:L
    lᵢ, lⱼ = l[i], l[j]
    ax = Makie.Axis(fig[i, j])
    Makie.lines!(ax, hs, getindex.(ts, i, j), color=color, linewidth=size, label="$lᵢ → $lⱼ")
    if i == j
      # show mean length
      Makie.lines!(ax, [zero(H), r[i]], [1.0, 0.0], color=color, linewidth=size, linestyle=:dash)
      # show proportion
      Makie.hlines!(ax, p[i], color=color, linewidth=size, linestyle=:dash)
    else
      # show base transiogram model
      Makie.lines!(ax, hs, getindex.(t̂s, i, j), color=color, linewidth=size, linestyle=:dot)
    end
    Makie.axislegend(position=i == j ? :rt : :rb)
  end
  fig
end
