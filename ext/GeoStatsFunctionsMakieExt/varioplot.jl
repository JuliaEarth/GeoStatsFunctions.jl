# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.@recipe(VarioPlot, γ) do scene
  Makie.Attributes(
    # common variogram options
    color=:slategray,
    size=1.5,
    maxlag=nothing,

    # empirical variogram options
    pointsize=12,
    showtext=true,
    textsize=12,
    showhist=true,
    histcolor=:slategray
  )
end

Makie.plottype(::EmpiricalVariogram) = VarioPlot{<:Tuple{EmpiricalVariogram}}

function Makie.plot!(plot::VarioPlot{<:Tuple{EmpiricalVariogram}})
  # retrieve variogram object
  γ = plot[:γ]

  # get the data
  x = Makie.@lift $γ.abscissas
  y = Makie.@lift $γ.ordinates
  n = Makie.@lift $γ.counts

  # discard empty bins
  x = Makie.@lift $x[$n .> 0]
  y = Makie.@lift $y[$n .> 0]
  n = Makie.@lift $n[$n .> 0]

  # visualize frequencies as bars
  if plot[:showhist][]
    f = Makie.@lift $n * (maximum($y) / maximum($n)) / 10
    Makie.barplot!(plot, x, f, color=plot[:histcolor], alpha=0.3, gap=0.0)
  end

  # visualize variogram
  Makie.scatterlines!(plot, x, y, color=plot[:color], markersize=plot[:pointsize], linewidth=plot[:size])

  # visualize text counts
  if plot[:showtext][]
    text = Makie.@lift string.($n)
    Makie.text!(plot, x, y, text=text, fontsize=plot[:textsize])
  end
end

Makie.plottype(::Variogram) = VarioPlot{<:Tuple{Variogram}}

function Makie.plot!(plot::VarioPlot{<:Tuple{Variogram}})
  # retrieve variogram object
  γ = plot[:γ]

  # retrieve maximum lag
  maxlag = plot[:maxlag]

  # adjust maximum lag
  H = if isnothing(maxlag[])
    Makie.@lift _maxlag($γ)
  else
    Makie.@lift _addunit($maxlag, u"m")
  end

  # lag range starting at 1e-6 to avoid nugget artifact
  hs = Makie.@lift range(1e-6unit($H), stop=$H, length=100)

  b = Makie.@lift metricball($γ)
  R = Makie.@lift rotation($b)
  r = Makie.@lift radii($b)
  n = Makie.@lift length($r)
  U = Makie.@lift eltype($r)
  p = Makie.@lift Point(ntuple(i -> $U(0), $n))
  for j in 1:n[]
    sⱼ = Dict(1 => :solid, 2 => :dash, 3 => :dot)[j]
    vⱼ = Makie.@lift $R * Vec(ntuple(i -> $U(i == j), $n))
    gs = Makie.@lift [$γ($p, $p + ustrip(h) * $vⱼ) for h in $hs]
    Makie.lines!(plot, hs, gs, color=plot[:color], linewidth=plot[:size], linestyle=sⱼ)
  end
end
