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

formula(::CarleEstimator, z₁ᵢ, z₁ⱼ, z₂ᵢ, z₂ⱼ) = SVector{2,Int}(z₁ᵢ * z₂ⱼ, z₁ᵢ)

normsum(::CarleEstimator, Σy, n) = (Σy[1] / Σy[2]) / n

combine(::CarleEstimator, yα, nα, yβ, nβ) = (yα * nα + yβ * nβ) / (nα + nβ)
