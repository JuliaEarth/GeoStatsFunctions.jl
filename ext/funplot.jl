# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function funplot(f; kwargs...)
  fig = Makie.Figure()
  funplot(fig, f; kwargs...)
end

function funplot(pos, f; link=nothing, kwargs...)
  # decide whether to link y axes of subplots
  ylink = isnothing(link) ? _ylink(f) : link

  # initialize layout
  n = nvariables(f)
  v = variables(f)
  layout = _layout(pos)
  axs = Makie.Axis[]
  for i in 1:n, j in 1:n
    issymmetric(f) && i < j && continue
    ax = Makie.Axis(layout[i, j])
    ax.title = issymmetric(f) ? "$(v[j]) → $(v[i])" : "$(v[i]) → $(v[j])"
    i == n && (ax.xlabel = "lag distance [m]")
    j == 1 && (ax.ylabel = _ylabel(f))
    i < n && Makie.hidexdecorations!(ax, grid=false)
    j > 1 && ylink && Makie.hideydecorations!(ax, grid=false)
    push!(axs, ax)
  end
  Makie.linkxaxes!(axs...)
  ylink && Makie.linkyaxes!(axs...)

  # fill layout with plots
  funplot!(layout, f; kwargs...)

  pos
end

function funplot!(
  layout::Makie.GridLayout,
  f::GeoStatsFunction;
  # common options
  color=:teal,
  linewidth=1.5,
  maxlag=nothing
)
  # maximum lag
  hmax = isnothing(maxlag) ? _maxlag(f) : GeoStatsFunctions.aslen(maxlag)

  # lag range starting at 1e-6 to avoid nugget artifact
  hs = range(1e-6 * oneunit(hmax), stop=hmax, length=100)

  # evaluate function
  F = _eval(f, hs)

  # aesthetic options
  label = Dict(1 => "1ˢᵗ axis", 2 => "2ⁿᵈ axis", 3 => "3ʳᵈ axis")
  linestyle = Dict(1 => :solid, 2 => :dash, 3 => :dot)

  # add plots to axes
  d = length(F)
  n = nvariables(f)
  for i in 1:n, j in 1:n
    issymmetric(f) && i < j && continue
    ax = Makie.content(layout[i, j])
    for (k, Fₖ) in enumerate(F)
      Makie.lines!(ax, ustrip.(u"m", hs), Fₖ[i, j]; color, linewidth, linestyle=linestyle[k], label=label[k])
    end
    position = (i == j) && isbanded(f) ? :rt : :rb
    d > 1 && Makie.axislegend(ax, position=position)
  end

  layout
end

function funplot!(
  layout::Makie.GridLayout,
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
  histcolor=color
)
  # maximum lag
  hmax = isnothing(maxlag) ? _maxlag(f) : GeoStatsFunctions.aslen(maxlag)

  # add plots to axes
  n = nvariables(f)
  for i in 1:n, j in 1:n
    issymmetric(f) && i < j && continue
    ax = Makie.content(layout[i, j])

    # retrieve coordinates and counts
    hs = f.abscissas
    fs = f.ordinates[i, j]
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

  layout
end

funplot!(fig::Makie.Figure, f; kwargs...) = (funplot!(fig.layout, f; kwargs...); fig)
