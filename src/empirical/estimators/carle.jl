# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CarleEstimator()

Carle's transiogram estimator (equation 10 of Carle, S.F. & Fogg, G.E. 1996).

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)
"""
struct CarleEstimator <: Estimator end
