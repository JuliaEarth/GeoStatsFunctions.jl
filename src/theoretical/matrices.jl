# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    baseratematrix(lengths, proportions)

Base rate matrix with given `lenghts` and `proportions`.

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)

* Carle et al 1998. [Conditional Simulation of Hydrofacies Architecture:
  A Transition Probability/Markov Approach](https://doi.org/10.2110/sepmcheg.01.147)
"""
function baseratematrix(l, p)
  nₗ = length(l)
  nₚ = length(p)

  # add units if necessary
  lᵤ = aslen.(l)

  # sanity checks
  if nₗ != nₚ
    throw(ArgumentError("lengths and proportions must have the same length"))
  end
  if !all(pᵢ -> 0 ≤ pᵢ ≤ 1, p)
    throw(ArgumentError("proportions must be in interval [0, 1]"))
  end
  if !(sum(p) ≈ 1)
    throw(ArgumentError("proportions must add up to unit"))
  end

  # Eq. 17 and 18 of Carle et al 1998.
  SMatrix{nₗ,nₗ}(i == j ? -1 / lᵤ[i] : (p[j] / (1 - p[i])) / lᵤ[i] for i in 1:nₗ, j in 1:nₗ)
end
