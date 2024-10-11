# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function transioplot(
  t::Transiogram;
  # common transiogram options
  color=:slategray,
  size=1.5,
  maxlag=nothing
)
  # effective ranges and labels
  r = GeoStatsFunctions.ranges(t)
  l = GeoStatsFunctions.levels(t)

  # number of labels
  L = length(l)

  # retrieve maximum lag
  H = if isnothing(maxlag)
    _maxlag(t)
  else
    _addunit(maxlag, u"m")
  end

  # transiogram up to maximum lag
  hs = range(zero(H), stop=H, length=100)
  ts = t.(hs)

  # relative proportions at h → ∞
  p = normalize(diag(t(100H)), 1)

  # base transiogram model
  t̂ = ExponentialTransiogram(r, p)
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
