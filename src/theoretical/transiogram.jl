# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Transiogram

A theoretical transiogram function (e.g. exponential transiogram).

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)
"""
abstract type Transiogram <: GeoStatsFunction end

# ---------------------
# GEOSTATSFUNCTION API
# ---------------------

isstationary(::Type{<:Transiogram}) = true

issymmetric(::Type{<:Transiogram}) = false

isbanded(::Type{<:Transiogram}) = true

function scale(t::Transiogram, s::Real)
  T = constructor(t)
  T(s * metricball(t); proportions=proportions(t))
end

nvariates(t::Transiogram) = length(proportions(t))

# ----------------
# TRANSIOGRAM API
# ----------------

meanlengths(t::Transiogram) = ntuple(i -> range(t), nvariates(t))

proportions(t::Transiogram) = t.proportions

# ----------------
# IMPLEMENTATIONS
# ----------------

include("transiogram/linear.jl")
include("transiogram/gaussian.jl")
include("transiogram/spherical.jl")
include("transiogram/exponential.jl")
include("transiogram/matrixexponential.jl")
include("transiogram/piecewiselinear.jl")
include("transiogram/carle.jl")
