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
    histcolor=:slategray,

    # empirical varioplane options
    colormap=:viridis,
    showrange=true,
    rangecolor=:white,
    rangemodel=SphericalVariogram
  )
end

# ----------
# EMPIRICAL
# ----------

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

Makie.plottype(::EmpiricalVarioplane) = VarioPlot{<:Tuple{EmpiricalVarioplane}}

function Makie.plot!(plot::VarioPlot{<:Tuple{EmpiricalVarioplane}})
  # retrieve varioplane object
  v = plot[:γ]

  # retrieve range model
  showrange = plot[:showrange]
  rangemodel = plot[:rangemodel]

  # underyling variograms
  γs = Makie.@lift $v.γs

  # polar angle
  θs = Makie.@lift $v.θs

  # polar radius
  rs = Makie.@lift ustrip.($γs[1].abscissas)

  # variogram ordinates for all variograms
  H = Makie.@lift let
    hs = map($γs) do γ
      # retrieve ordinates without units
      ys = ustrip.(γ.ordinates)

      # handle NaN ordinates (i.e. empty bins)
      isnan(ys[1]) && (ys[1] = 0)
      for i in 2:length(ys)
        isnan(ys[i]) && (ys[i] = ys[i - 1])
      end

      ys
    end
    reduce(hcat, hs)
  end

  # exploit symmetry
  θs = Makie.@lift range(0, 2π, length=2 * length($θs))
  H = Makie.@lift [$H $H]

  # hide hole at center
  rs = Makie.@lift [0; $rs]
  H = Makie.@lift [$H[1:1, :]; $H]

  # transpose for plotting
  H = Makie.@lift transpose($H)

  Makie.surface!(plot, θs, rs, H, colormap=plot[:colormap], shading=Makie.NoShading)

  # show model range
  if showrange[]
    ls = Makie.@lift [ustrip(range(GeoStatsFunctions.fit($rangemodel, γ))) for γ in $γs]
    ls = Makie.@lift [$ls; $ls]
    hs = Makie.@lift fill(maximum($H) + 1, length($ls))
    Makie.lines!(plot, θs, ls, hs, color=plot[:rangecolor])
  end
end

# ------------
# THEORETICAL
# ------------

Makie.plottype(::Variogram) = VarioPlot{<:Tuple{Variogram}}

function Makie.plot!(plot::VarioPlot{<:Tuple{Variogram}})
  # retrieve variogram object
  γ = plot[:γ]

  # retrieve maximum lag
  maxlag = plot[:maxlag]

  H = if isnothing(maxlag[])
    Makie.@lift _maxlag($γ)
  else
    Makie.@lift _addunit($maxlag, u"m")
  end

  # start at 1e-6 instead of 0 to avoid
  # nugget artifact in visualization
  x = Makie.@lift range(1e-6unit($H), stop=$H, length=100)
  y = Makie.@lift $γ.($x)

  # visualize variogram
  Makie.lines!(plot, x, y, color=plot[:color], linewidth=plot[:size])
end
