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

nvariates(t::Transiogram) = length(proportions(t))

# ----------------
# TRANSIOGRAM API
# ----------------

meanlengths(t::Transiogram) = radii(metricball(t))

proportions(t::Transiogram) = Tuple(normalize(diag(t(100range(t))), 1))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("transiogram/linear.jl")
include("transiogram/gaussian.jl")
include("transiogram/spherical.jl")
include("transiogram/exponential.jl")
include("transiogram/matrixexponential.jl")
include("transiogram/piecewiselinear.jl")
