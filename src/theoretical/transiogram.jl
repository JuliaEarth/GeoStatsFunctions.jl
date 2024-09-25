# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Transiogram

A theoretical transiogram model (e.g. exponential transiogram).

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)
"""
abstract type Transiogram end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("transiogram/exponential.jl")
