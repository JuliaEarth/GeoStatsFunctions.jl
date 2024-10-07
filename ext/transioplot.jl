# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function transioplot(t::Transiogram;
  # common transiogram options
  color=:slategray,
  size=1.5,
  maxlag=nothing
)
  # retrieve maximum lag
  H = if isnothing(maxlag)
    _maxlag(t)
  else
    _addunit(maxlag, u"m")
  end

  # effective ranges
  r = GeoStatsFunctions.ranges(t)

  # transiogram up to maximum lag
  x = range(zero(H), stop=H, length=100)
  y = t.(x)

  # number of levels
  L = Base.size(first(y), 1)

  fig = Makie.Figure()
  for i in 1:L, j in 1:L
    ax = Makie.Axis(fig[i, j])
    ys = getindex.(y, i, j)
    Makie.lines!(ax, x, ys, color=color, linewidth=size, label = "$i → $j")
    if i == j
      # effective range
      Makie.lines!(ax, [zero(H), r[i]], [1.0, 0.0], color=color, linewidth=size, linestyle=:dash)
    end
    Makie.axislegend(position = i == j ? :rt : :rb)
  end
  fig
end
