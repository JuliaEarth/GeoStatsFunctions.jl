# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Variogram

A theoretical variogram function (e.g. Gaussian variogram).
"""
abstract type Variogram <: GeoStatsFunction end

"""
    isstationary(γ)

Check if variogram `γ` possesses the 2nd-order stationary property.
"""
isstationary(γ::Variogram) = isstationary(typeof(γ))

"""
    sill(γ)

Return the sill of the variogram `γ`.
"""
sill(γ::Variogram) = γ.sill

"""
    nugget(γ)

Return the nugget of the variogram `γ`.
"""
nugget(γ::Variogram) = γ.nugget

"""
    scale(γ, s)

Scale metric ball of variogram `γ` with strictly
positive scaling factor `s`.
"""
function scale(γ::Variogram, s::Real)
  V = constructor(γ)
  V(s * metricball(γ); sill=sill(γ), nugget=nugget(γ))
end

"""
    structures(γ)

Return the individual structures of a (possibly nested)
variogram as a tuple. The structures are the total nugget
`cₒ`, and the coefficients (or contributions) `c[i]` for the
remaining non-trivial structures `g[i]` after normalization
(i.e. sill=1, nugget=0).

## Examples

```julia
γ₁ = GaussianVariogram(nugget=1, sill=2)
γ₂ = SphericalVariogram(nugget=2, sill=3)

γ = 2γ₁ + 3γ₂

cₒ, c, g = structures(γ)
```
"""
function structures(γ::Variogram)
  cₒ = nugget(γ)
  c = sill(γ) - nugget(γ)
  T = typeof(c)
  γ = @set γ.sill = one(T)
  γ = @set γ.nugget = zero(T)
  cₒ, (c,), (γ,)
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, γ::T) where {T<:Variogram}
  name = string(nameof(T))
  _showcompact(io, name, γ)
end

function Base.show(io::IO, ::MIME"text/plain", γ::T) where {T<:Variogram}
  name = string(nameof(T))
  _showfull(io, name, γ)
end

function _showcompact(io, name, γ::T) where {T<:Variogram}
  params = String[]
  for fn in fieldnames(T)
    val = getfield(γ, fn)
    if val isa MetricBall
      if isisotropic(val)
        r = first(radii(val))
        push!(params, "range: $r")
      else
        r = Tuple(radii(val))
        push!(params, "ranges: $r")
      end
      d = nameof(typeof(metric(val)))
      push!(params, "distance: $d")
    else
      push!(params, "$fn: $val")
    end
  end
  print(io, name, "(", join(params, ", "), ")")
end

function _showfull(io, name, γ::T) where {T<:Variogram}
  header = isisotropic(γ) ? name : name * " (anisotropic)"
  params = String[]
  fnames = fieldnames(T)
  len = length(fnames)
  for (i, fn) in enumerate(fnames)
    div = i == len ? "└─" : "├─"
    val = getfield(γ, fn)
    if val isa MetricBall
      if isisotropic(val)
        r = first(radii(val))
        push!(params, "├─ range: $r")
      else
        r = Tuple(radii(val))
        push!(params, "└─ ranges: $r")
      end
      m = nameof(typeof(metric(val)))
      push!(params, "$div distance: $m")
    else
      push!(params, "$div $(fn): $val")
    end
  end
  println(io, header)
  print(io, join(params, "\n"))
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("variogram/gaussian.jl")
include("variogram/exponential.jl")
include("variogram/spherical.jl")
include("variogram/matern.jl")
include("variogram/cubic.jl")
include("variogram/pentaspherical.jl")
include("variogram/sinehole.jl")
include("variogram/power.jl")
include("variogram/nugget.jl")
include("variogram/circular.jl")
include("variogram/nested.jl")
