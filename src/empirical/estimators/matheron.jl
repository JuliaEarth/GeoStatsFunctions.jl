# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MatheronEstimator()

Matheron's variogram estimator.
"""
struct MatheronEstimator <: Estimator end

accumterm(::MatheronEstimator, z₁ᵢ, z₁ⱼ, z₂ᵢ, z₂ⱼ) = (z₁ᵢ - z₁ⱼ) ⋅ (z₂ᵢ - z₂ⱼ)

accumnorm(::MatheronEstimator, Σy, n) = Σy / (2n)

mergerule(::MatheronEstimator, yα, nα, yβ, nβ) = (yα * nα + yβ * nβ) / (nα + nβ)
