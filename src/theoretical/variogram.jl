# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Variogram

A theoretical variogram function (e.g. Gaussian variogram).
"""
abstract type Variogram <: GeoStatsFunction end

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

# ----------------
# IMPLEMENTATIONS
# ----------------

include("variogram/gaussian.jl")
include("variogram/spherical.jl")
include("variogram/exponential.jl")
include("variogram/matern.jl")
include("variogram/cubic.jl")
include("variogram/pentaspherical.jl")
include("variogram/sinehole.jl")
include("variogram/circular.jl")
include("variogram/power.jl")
include("variogram/nugget.jl")
include("variogram/nested.jl")
