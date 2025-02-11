# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function funplot(
  f::GeoStatsFunction;
  # common options
  color=:teal,
  size=1.5,
  maxlag=nothing,
  labels=nothing
)
  # auxiliary parameters
  n = nvariates(f)
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

  # maximum lag
  hmax = isnothing(maxlag) ? _maxlag(f) : _addunit(maxlag, u"m")

  # lag range starting at 1e-6 to avoid nugget artifact
  hs = range(1e-6 * oneunit(hmax), stop=hmax, length=100)

  # evaluate main function
  F = map(v) do vⱼ
    fs = [f(p, p + ustrip(h) * vⱼ) for h in hs]
    [getindex.(fs, i, j) for i in 1:n, j in 1:n]
  end

  # evaluate base function
  λ, π, F̂ = if _istransiogram(f)
    λ = meanlengths(f)
    π = proportions(f)
    f̂ = MatrixExponentialTransiogram(lengths=λ, proportions=π)
    F̂ = map(v) do vⱼ
      fs = [f̂(p, p + ustrip(h) * vⱼ) for h in hs]
      [getindex.(fs, i, j) for i in 1:n, j in 1:n]
    end
    λ, π, F̂
  else
    nothing, nothing, nothing
  end

  # aesthetic options
  label = Dict(1 => "1ˢᵗ axis", 2 => "2ⁿᵈ axis", 3 => "3ʳᵈ axis")
  style = Dict(1 => :solid, 2 => :dash, 3 => :dot)
  vars = isnothing(labels) ? (1:n) : labels

  # build figure
  fig = Makie.Figure()
  for i in 1:n, j in 1:n
    title = n > 1 ? "$(vars[i]) → $(vars[j])" : ""
    ax = Makie.Axis(fig[i, j], title=title)
    for (k, Fₖ) in enumerate(F)
      # display main function
      Makie.lines!(ax, hs, Fₖ[i, j], color=color, linewidth=size, linestyle=style[k], label=label[k])
      if _istransiogram(f)
        if i == j
          # display mean lengths
          Makie.lines!(ax, [zero(hmax), λ[i]], [1.0, 0.0], color=:slategray, linewidth=size, linestyle=:dash)
          # display proportions
          Makie.hlines!(ax, π[i], color=:slategray, linewidth=size, linestyle=:dash)
        else
          # display base function
          F̂ₖ = F̂[k]
          Makie.lines!(ax, hs, F̂ₖ[i, j], color=:slategray, linewidth=size, linestyle=:dot)
        end
      end
    end
    position = (i == j) && isbanded(f) ? :rt : :rb
    d > 1 && Makie.axislegend(position=position)
  end
  fig
end

function funplot(
  f::EmpiricalGeoStatsFunction;
  # common options
  color=:teal,
  size=1.5,
  maxlag=nothing,
  labels=nothing,
  # empirical options
  pointsize=12,
  showtext=true,
  textsize=12,
  showhist=true,
  histcolor=:teal
)
  # number of variables
  n = nvariates(f)

  # aesthetic options
  vars = isnothing(labels) ? (1:n) : labels

  fig = Makie.Figure()
  for i in 1:n, j in 1:n
    title = n > 1 ? "$(vars[i]) → $(vars[j])" : ""
    ax = Makie.Axis(fig[i, j], title=title)

    # retrieve coordinates and counts
    hs = f.abscissas
    fs = _istransiogram(f) ? f.ordinates[i, j] : f.ordinates
    ns = f.counts

    # retrieve maximum lag
    hmax = isnothing(maxlag) ? _maxlag(f) : _addunit(maxlag, u"m")

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
      Makie.barplot!(ax, hs, ms, color=histcolor, alpha=0.3, gap=0.0)
    end

    # visualize function
    Makie.scatterlines!(ax, hs, fs, color=color, markersize=pointsize, linewidth=size)

    # visualize text counts
    if showtext
      Makie.text!(ax, hs, fs, text=string.(ns), fontsize=textsize)
    end
  end
  fig
end
