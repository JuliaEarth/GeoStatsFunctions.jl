# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

funplot(f; kwargs...) = _funplot(f; kwargs...)

function _funplot(
  f::GeoStatsFunction;

  # common options
  color=:slategray,
  size=1.5,
  maxlag=nothing,
  levels=nothing
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

  # evaluate function
  F = map(v) do vⱼ
    fs = [f(p, p + ustrip(h) * vⱼ) for h in hs]
    [getindex.(fs, i, j) for i in 1:n, j in 1:n]
  end

  # aesthetic options
  label = Dict(1 => "1ˢᵗ axis", 2 => "2ⁿᵈ axis", 3 => "3ʳᵈ axis")
  style = Dict(1 => :solid, 2 => :dash, 3 => :dot)
  level = isnothing(levels) ? (1:n) : levels

  # build figure
  fig = Makie.Figure()
  for i in 1:n, j in 1:n
    lᵢ, lⱼ = level[i], level[j]
    title = n > 1 ? "$lᵢ → $lⱼ" : ""
    ax = Makie.Axis(fig[i, j], title=title)
    for (k, Fₖ) in enumerate(F)
      Makie.lines!(ax, hs, Fₖ[i, j], color=color, linewidth=size, linestyle=style[k], label=label[k])
    end
    position = (i == j) && isbanded(f) ? :rt : :rb
    d > 1 && Makie.axislegend(position=position)
  end
  fig
end

# ----------------
# SPECIALIZATIONS
# ----------------

include("funplot/variogram.jl")
include("funplot/transiogram.jl")
