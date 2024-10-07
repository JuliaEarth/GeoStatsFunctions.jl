# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.@recipe(VarioPlot, γ) do scene
  Makie.Attributes(
    # empirical variogram options
    color=:slategray,
    pointsize=12,
    segmentsize=1.5,
    showtext=true,
    textsize=12,
    showhist=true,
    histcolor=:slategray,

    # empirical varioplane options
    colormap=:viridis,
    showrange=true,
    rangecolor=:white,
    rangemodel=SphericalVariogram,

    # theoretical variogram options
    maxlag=nothing
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
  xyn = Makie.@lift values($γ)
  x = Makie.@lift $xyn[1]
  y = Makie.@lift $xyn[2]
  n = Makie.@lift $xyn[3]

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
  Makie.scatterlines!(plot, x, y, color=plot[:color], markersize=plot[:pointsize], linewidth=plot[:segmentsize])

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
  rs = Makie.@lift ustrip.(values($γs[1])[1])

  # variogram values for all variograms
  Z = Makie.@lift let
    zs = map($γs) do γ
      zs = ustrip.(values(γ)[2])

      # handle NaN values (i.e. empty bins)
      isnan(zs[1]) && (zs[1] = 0)
      for i in 2:length(zs)
        isnan(zs[i]) && (zs[i] = zs[i - 1])
      end

      zs
    end
    reduce(hcat, zs)
  end

  # exploit symmetry
  θs = Makie.@lift range(0, 2π, length=2 * length($θs))
  Z = Makie.@lift [$Z $Z]

  # hide hole at center
  rs = Makie.@lift [0; $rs]
  Z = Makie.@lift [$Z[1:1, :]; $Z]

  # transpose for plotting
  Z = Makie.@lift transpose($Z)

  Makie.surface!(plot, θs, rs, Z, colormap=plot[:colormap], shading=Makie.NoShading)

  # show model range
  if showrange[]
    ls = Makie.@lift [ustrip(range(GeoStatsFunctions.fit($rangemodel, γ))) for γ in $γs]
    ls = Makie.@lift [$ls; $ls]
    zs = Makie.@lift fill(maximum($Z) + 1, length($ls))
    Makie.lines!(plot, θs, ls, zs, color=plot[:rangecolor])
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
  Makie.lines!(plot, x, y, color=plot[:color])
end

# ---------------------------
# helper types and functions
# ---------------------------

const Len{T} = Quantity{T,u"𝐋"}

_maxlag(γ::Variogram) = 3range(γ)
_maxlag(::PowerVariogram) = 3.0u"m"
_maxlag(::NuggetEffect) = 3.0u"m"

_addunit(x::Number, u) = x * u
_addunit(x::Len, _) = x
_addunit(x::Quantity, _) = throw(ArgumentError("$(unit(x)) is not a valid length unit"))
