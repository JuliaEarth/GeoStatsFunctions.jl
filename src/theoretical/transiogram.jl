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

Base.range(t::Transiogram) = maximum(ranges(t))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("transiogram/exponential.jl")
include("transiogram/piecewiselinear.jl")
