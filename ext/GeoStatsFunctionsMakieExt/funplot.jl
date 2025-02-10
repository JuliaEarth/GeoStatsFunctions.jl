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
  # evaluate at rotated basis
  hs, N, F = _evalbasis(f, maxlag)

  # aesthetic options
  label = Dict(1 => "1st axis", 2 => "2nd axis", 3 => "3rd axis")
  style = Dict(1 => :solid, 2 => :dash, 3 => :dot)
  level = isnothing(levels) ? (1:N) : levels

  fig = Makie.Figure()
  for i in 1:N, j in 1:N
    lᵢ, lⱼ = level[i], level[j]
    title = N > 1 ? "$lᵢ → $lⱼ" : ""
    ax = Makie.Axis(fig[i, j], title = title)
    for (k, Fₖ) in enumerate(F)
      Makie.lines!(ax, hs, Fₖ[i, j], color=color, linewidth=size, linestyle=style[k], label=label[k])
    end
    position = (i == j) && isbanded(f) ? :rt : :rb
    length(F) > 1 && Makie.axislegend(position=position)
  end
  fig
end

# ----------------
# SPECIALIZATIONS
# ----------------

include("funplot/variogram.jl")
include("funplot/transiogram.jl")

# -----------------
# HELPER FUNCTIONS
# -----------------

function _evalbasis(f, maxlag)
  # auxiliary parameters
  b = metricball(f)
  R = rotation(b)
  r = radii(b)
  n = length(r)
  U = eltype(r)

  # reference point
  p = Point(ntuple(i -> U(0), n))

  # rotated basis
  v = ntuple(n) do j
    R * Vec(ntuple(i -> U(i == j), n))
  end

  # maximum lag
  H = isnothing(maxlag) ? _maxlag(f) : _addunit(maxlag, u"m")

  # lag range starting at 1e-6 to avoid nugget artifact
  hs = range(1e-6 * oneunit(H), stop=H, length=100)

  # number of variates
  N = nvariates(f)

  # reshape results
  F = map(v) do vⱼ
    fs = [f(p, p + ustrip(h) * vⱼ) for h in hs]
    [getindex.(fs, i, j) for i in 1:N, j in 1:N]
  end

  hs, N, F
end
