# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function transioplot(t::Transiogram;
  # empirical transiogram options
  color=:slategray,
  pointsize=12,
  segmentsize=1.5,
  showtext=true,
  textsize=12,
  showhist=true,
  histcolor=:slategray,

  # theoretical transiogram options
  maxlag=nothing
)
  # retrieve maximum lag
  H = if isnothing(maxlag)
    _maxlag(t)
  else
    _addunit(maxlag, u"m")
  end

  # transiogram up to maximum lag
  x = range(zero(H), stop=H, length=100)
  y = t.(x)

  # number of levels
  L = size(first(y), 1)

  fig = Makie.Figure()
  for i in 1:L, j in 1:L
    ax = Makie.Axis(fig[i, j])
    ys = getindex.(y, i, j)
    Makie.lines!(ax, x, ys, color = color)
  end
  fig
end
