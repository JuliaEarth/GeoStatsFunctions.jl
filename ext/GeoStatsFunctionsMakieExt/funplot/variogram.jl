# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _funplot(
  γ::EmpiricalVariogram;

  # common options
  color=:slategray,
  size=1.5,
  maxlag=nothing,

  # empirical options
  pointsize=12,
  showtext=true,
  textsize=12,
  showhist=true,
  histcolor=:slategray
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

function _funplot(
  γ::Variogram;

  # common options
  color=:slategray,
  size=1.5,
  maxlag=nothing
)
  # auxiliary parameters
  b = metricball(γ)
  R = rotation(b)
  r = radii(b)
  n = length(r)
  U = eltype(r)

  # reference point
  p = Point(ntuple(i -> U(0), n))

  # retrieve maximum lag
  H = isnothing(maxlag) ? _maxlag(γ) : _addunit(maxlag, u"m")

  # lag range starting at 1e-6 to avoid nugget artifact
  hs = range(1e-6 * oneunit(H), stop=H, length=100)

  # initialize figure and axis
  fig = Makie.Figure()
  ax = Makie.Axis(fig[1, 1])
  for j in 1:n
    sⱼ = Dict(1 => :solid, 2 => :dash, 3 => :dot)[j]
    lⱼ = Dict(1 => "1st axis", 2 => "2nd axis", 3 => "3rd axis")[j]
    vⱼ = R * Vec(ntuple(i -> U(i == j), n))
    gs = [γ(p, p + ustrip(h) * vⱼ) for h in hs]
    Makie.lines!(ax, hs, gs, color=color, linewidth=size, linestyle=sⱼ, label=lⱼ)
  end
  n > 1 && Makie.axislegend(position=:rb)

  fig
end
