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
  # effective ranges and variable levels
  r = GeoStatsFunctions.ranges(t)
  l = GeoStatsFunctions.levels(t)

  # number of levels
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

  # approximate transiogram at h → ∞
  t∞ = t(100H)

  fig = Makie.Figure()
  for i in 1:L, j in 1:L
    lᵢ, lⱼ = l[i], l[j]
    ax = Makie.Axis(fig[i, j])
    ys = getindex.(ts, i, j)
    Makie.lines!(ax, hs, ys, color=color, linewidth=size, label="$lᵢ → $lⱼ")
    if i == j
      # show mean length
      Makie.lines!(ax, [zero(H), r[i]], [1.0, 0.0], color=color, linewidth=size, linestyle=:dash)
      # show proportion
      Makie.hlines!(ax, t∞[i, j], color=color, linewidth=size, linestyle=:dash)
    else
      # show random transition model
      pᵢ = t∞[i, i]
      pⱼ = t∞[j, j]
      tᵢᵢ = getindex.(ts, i, i)
      rs = @. (1 - tᵢᵢ) * pⱼ / (1 - pᵢ)
      Makie.lines!(ax, hs, rs, color=color, linewidth=size, linestyle=:dot)
    end
    Makie.axislegend(position=i == j ? :rt : :rb)
  end
  fig
end
