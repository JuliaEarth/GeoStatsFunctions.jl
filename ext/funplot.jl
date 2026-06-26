# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function funplot(f; names=nothing, kwargs...)
  # initialize figure
  n = nvariates(f)
  v = isnothing(names) ? (1:n) : names
  fig = Makie.Figure()
  for i in 1:n, j in 1:n
    ax = Makie.Axis(fig[i, j])
    ax.title = n > 1 ? "$(v[i]) → $(v[j])" : ""
    ax.xlabel = i == n ? "lag distance [m]" : ""
  end
  Makie.linkaxes!(fig.content...)

  # fill figure with plots
  funplot!(fig, f; kwargs...)
end

function funplot!(
  fig::Makie.Figure,
  f::GeoStatsFunction;
  # common options
  color=:teal,
  linewidth=1.5,
  maxlag=nothing,
  # theoretical options
  showmeanlengths=false,
  showproportions=true
)
  # maximum lag
  hmax = isnothing(maxlag) ? _maxlag(f) : GeoStatsFunctions.aslen(maxlag)

  # lag range starting at 1e-6 to avoid nugget artifact
  hs = range(1e-6 * oneunit(hmax), stop=hmax, length=100)

  # evaluate function
  F = _eval(f, hs)

  # mean lengths and proportions
  l = _istransiogram(f) ? meanlengths(f) : nothing
  p = _istransiogram(f) ? proportions(f) : nothing

  # aesthetic options
  label = Dict(1 => "1ˢᵗ axis", 2 => "2ⁿᵈ axis", 3 => "3ʳᵈ axis")
  linestyle = Dict(1 => :solid, 2 => :dash, 3 => :dot)

  # add plots to axes
  d = length(F)
  n = nvariates(f)
  I = LinearIndices((n, n))
  for i in 1:n, j in 1:n
    ax = fig.content[I[j, i]]
    for (k, Fₖ) in enumerate(F)
      # display function
      Makie.lines!(ax, ustrip.(u"m", hs), Fₖ[i, j]; color, linewidth, linestyle=linestyle[k], label=label[k])
      if _istransiogram(f)
        if i == j
          if showmeanlengths
            # display mean lengths
            Makie.lines!(
              ax,
              ustrip.(u"m", [zero(hmax), l[i]]),
              [1.0, 0.0];
              color=:slategray,
              linewidth,
              linestyle=:dash
            )
          end
          if showproportions
            # display proportions
            Makie.hlines!(ax, p[i]; color=:slategray, linewidth, linestyle=:dash)
          end
        end
      end
    end
    position = (i == j) && isbanded(f) ? :rt : :rb
    d > 1 && Makie.axislegend(ax, position=position)
  end

  fig
end

_eval(f, hs) = isisotropic(f) ? _isoeval(f, hs) : _anisoeval(f, hs)

function _isoeval(f, hs)
  # auxiliary parameters
  n = nvariates(f)

  # evaluate all lags
  fs = f.(hs)

  # reshape result
  Fᵢₛₒ = [getindex.(fs, i, j) for i in 1:n, j in 1:n]

  [Fᵢₛₒ]
end

function _anisoeval(f, hs)
  # auxiliary parameters
  n = nvariates(f)

  # reference point and basis vectors
  p, v = _anisobasis(f)

  # evaluate along basis vectors
  map(v) do vⱼ
    fs = [f(p, p + ustrip(h) * vⱼ) for h in hs]
    [getindex.(fs, i, j) for i in 1:n, j in 1:n]
  end
end

function _anisobasis(f)
  # auxiliary parameters
  b = metricball(f)
  R = rotation(b)
  r = radii(b)
  d = length(r)
  U = eltype(r)

  # reference point
  p = Point(ntuple(i -> U(0), d))

  # rotated basis
  v = ntuple(d) do j
    R * Vec(ntuple(i -> U(i == j), d))
  end

  p, v
end

function _anisobasis(f::CarleTransiogram{N}) where {N}
  # auxiliary parameters
  U = typeof(range(f))

  # reference point
  p = Point(ntuple(i -> U(0), N))

  # axis-aligned basis
  v = ntuple(N) do j
    Vec(ntuple(i -> U(i == j), N))
  end

  p, v
end

function funplot!(
  fig::Makie.Figure,
  f::EmpiricalGeoStatsFunction;
  # common options
  color=:slategray,
  linewidth=1.5,
  maxlag=nothing,
  # empirical options
  linestyle=:solid,
  pointsize=12,
  showtext=true,
  textsize=12,
  showhist=true,
  histcolor=:slategray
)
  # maximum lag
  hmax = isnothing(maxlag) ? _maxlag(f) : GeoStatsFunctions.aslen(maxlag)

  # add plots to axes
  n = nvariates(f)
  I = LinearIndices((n, n))
  for i in 1:n, j in 1:n
    ax = fig.content[I[j, i]]

    # retrieve coordinates and counts
    hs = f.abscissas
    fs = _istransiogram(f) ? f.ordinates[i, j] : f.ordinates
    ns = f.counts

    # discard empty bins
    hs = hs[ns .> 0]
    fs = fs[ns .> 0]
    ns = ns[ns .> 0]

    # discard above maximum lag
    ns = ns[hs .≤ hmax]
    fs = fs[hs .≤ hmax]
    hs = hs[hs .≤ hmax]

    # visualize histogram
    if showhist
      ms = ns * (maximum(fs) / maximum(ns)) / 10
      Makie.barplot!(ax, ustrip.(u"m", hs), ms; color=histcolor, alpha=0.3, gap=0.0)
    end

    # visualize function
    Makie.scatterlines!(ax, ustrip.(u"m", hs), fs; color, markersize=pointsize, linewidth, linestyle)

    # visualize text counts
    if showtext
      Makie.annotation!(ax, ustrip.(u"m", hs), fs; text=string.(ns), fontsize=textsize)
    end
  end

  fig
end
