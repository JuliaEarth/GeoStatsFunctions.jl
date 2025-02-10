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

  # transiogram options
  levels=nothing
)
  # setup rotated basis
  p, v, hs = _setupbasis(f, maxlag)

  # aesthetic options
  label = Dict(1 => "1st axis", 2 => "2nd axis", 3 => "3rd axis")
  style = Dict(1 => :solid, 2 => :dash, 3 => :dot)
  aes = (; color, size, label, style, levels)

  # build figure
  _buildfigure(f, p, v, hs, aes)
end

# ----------------
# SPECIALIZATIONS
# ----------------

include("funplot/variogram.jl")
include("funplot/covariance.jl")
include("funplot/transiogram.jl")

# -----------------
# HELPER FUNCTIONS
# -----------------

function _setupbasis(f, maxlag)
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

  p, v, hs
end

function _buildfigure(f, p, v, hs, aes)
  fig = Makie.Figure()
  ax = Makie.Axis(fig[1, 1])
  for (j, vⱼ) in enumerate(v)
    fs = [f(p, p + ustrip(h) * vⱼ) for h in hs]
    Makie.lines!(ax, hs, fs, color=aes.color, linewidth=aes.size, linestyle=aes.style[j], label=aes.label[j])
  end
  length(v) > 1 && Makie.axislegend(position=:rb)
  fig
end
