# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function funplot(f; labels=nothing, kwargs...)
  fig = Makie.Figure()
  n = nvariates(f)
  l = isnothing(labels) ? (1:n) : labels
  for i in 1:n, j in 1:n
    title = n > 1 ? "$(l[i]) → $(l[j])" : ""
    Makie.Axis(fig[i, j], title=title)
  end
  funplot!(fig, f; kwargs...)
  fig
end

function funplot(gp::GridPos, f; labels=nothing, kwargs...)
  n = nvariates(f)
  l = isnothing(labels) ? (1:n) : labels
  layout = _layout(gp)
  for i in 1:n, j in 1:n
    title = n > 1 ? "$(l[i]) → $(l[j])" : ""
    Makie.Axis(layout[i, j], title=title)
  end
  funplot!(layout, f; kwargs...)
  gp
end

function funplot!(
  layout::Union{Makie.Figure, Makie.GridLayout},
  f::GeoStatsFunction;
  color=:teal,
  size=1.5,
  maxlag=nothing
)
  hmax = isnothing(maxlag) ? _maxlag(f) : GeoStatsFunctions.aslen(maxlag)
  hs = range(1e-6 * oneunit(hmax), stop=hmax, length=100)
  F = _eval(f, hs)
  l = _istransiogram(f) ? meanlengths(f) : nothing
  p = _istransiogram(f) ? proportions(f) : nothing
  label = Dict(1 => "1ˢᵗ axis", 2 => "2ⁿᵈ axis", 3 => "3ʳᵈ axis")
  style = Dict(1 => :solid, 2 => :dash, 3 => :dot)
  d = length(F)
  n = nvariates(f)
  gl = layout isa Makie.Figure ? layout.layout : layout
  I = LinearIndices((n, n))
  for i in 1:n, j in 1:n
    ax = gl.content[I[j, i]].content
    for (k, Fₖ) in enumerate(F)
      Makie.lines!(ax, hs, Fₖ[i, j], color=color, linewidth=size, linestyle=style[k], label=label[k])
      if _istransiogram(f)
        if i == j
          Makie.lines!(ax, [zero(hmax), l[i]], [1.0, 0.0], color=:slategray, linewidth=size, linestyle=:dash)
          Makie.hlines!(ax, p[i], color=:slategray, linewidth=size, linestyle=:dash)
        end
      end
    end
    position = (i == j) && isbanded(f) ? :rt : :rb
    d > 1 && Makie.axislegend(ax, position=position)
  end
  layout
end

_eval(f, hs) = isisotropic(f) ? _isoeval(f, hs) : _anisoeval(f, hs)

function _isoeval(f, hs)
  n = nvariates(f)
  fs = f.(hs)
  Fᵢₛₒ = [getindex.(fs, i, j) for i in 1:n, j in 1:n]
  [Fᵢₛₒ]
end

function _anisoeval(f, hs)
  n = nvariates(f)
  p, v = _anisobasis(f)
  map(v) do vⱼ
    fs = [f(p, p + ustrip(h) * vⱼ) for h in hs]
    [getindex.(fs, i, j) for i in 1:n, j in 1:n]
  end
end

function _anisobasis(f)
  b = metricball(f)
  R = rotation(b)
  r = radii(b)
  d = length(r)
  U = eltype(r)
  p = Point(ntuple(i -> U(0), d))
  v = ntuple(d) do j
    R * Vec(ntuple(i -> U(i == j), d))
  end
  p, v
end

function _anisobasis(f::CarleTransiogram{N}) where {N}
  U = typeof(range(f))
  p = Point(ntuple(i -> U(0), N))
  v = ntuple(N) do j
    Vec(ntuple(i -> U(i == j), N))
  end
  p, v
end

function funplot!(
  layout::Union{Makie.Figure, Makie.GridLayout},
  f::EmpiricalGeoStatsFunction;
  color=:slategray,
  size=1.5,
  maxlag=nothing,
  style=:solid,
  pointsize=12,
  showtext=true,
  textsize=12,
  showhist=true,
  histcolor=:slategray
)
  n = nvariates(f)
  gl = layout isa Makie.Figure ? layout.layout : layout
  I = LinearIndices((n, n))
  for i in 1:n, j in 1:n
    ax = gl.content[I[j, i]].content
    hs = f.abscissas
    fs = _istransiogram(f) ? f.ordinates[i, j] : f.ordinates
    ns = f.counts
    hmax = isnothing(maxlag) ? _maxlag(f) : GeoStatsFunctions.aslen(maxlag)
    hs = hs[ns .> 0]; fs = fs[ns .> 0]; ns = ns[ns .> 0]
    ns = ns[hs .≤ hmax]; fs = fs[hs .≤ hmax]; hs = hs[hs .≤ hmax]
    if showhist
      ms = ns * (maximum(fs) / maximum(ns)) / 10
      Makie.barplot!(ax, hs, ms, color=histcolor, alpha=0.3, gap=0.0)
    end
    Makie.scatterlines!(ax, hs, fs, color=color, markersize=pointsize, linewidth=size, linestyle=style)
    if showtext
      Makie.annotation!(ax, hs, fs, text=string.(ns), fontsize=textsize)
    end
  end
  layout
end