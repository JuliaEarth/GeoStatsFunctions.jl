# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ----------
# EMPIRICAL
# ----------

"""
    ratematrix(data, var; [parameters])

Estimate the transition rate matrix for categorical variable
`var` stored in geospatial `data` using a minimum lag `minlag`.

See also [`probmatrix`](@ref), [`countmatrix`](@ref)
"""
ratematrix(data::AbstractGeoTable, var; minlag=_defaultminlag(data)) = _ratematrix(_domvals(data, var)..., minlag)

"""
    probmatrix(data, var; [parameters])

Estimate the transition probability matrix for categorical variable
`var` stored in geospatial `data` using a minimum lag `minlag`.

See also [`ratematrix`](@ref), [`countmatrix`](@ref)
"""
probmatrix(data::AbstractGeoTable, var; minlag=_defaultminlag(data)) = _probmatrix(_domvals(data, var)..., minlag)

"""
    countmatrix(data, var; [parameters])

Estimate the transition count matrix for categorical variable
`var` stored in geospatial `data` using a minimum lag `minlag`.

See also [`ratematrix`](@ref), [`probmatrix`](@ref)
"""
countmatrix(data::AbstractGeoTable, var; minlag=_defaultminlag(data)) = _countmatrix(_domvals(data, var)..., minlag)

function _ratematrix(dom, vals, minlag)
  # transition probabilities
  T = _probmatrix(dom, vals, minlag)

  # convert probabilities into rates
  δh = minlag
  log(T) / δh
end

function _probmatrix(dom, vals, minlag)
  # transition counts
  C = _countmatrix(dom, vals, minlag)

  # convert counts into probabilities
  C ./ sum(C, dims=2)
end

function _countmatrix(dom, vals, minlag)
  # sanity check
  if elscitype(vals) != Categorical
    throw(ArgumentError("count matrix only defined for categorical variables"))
  end

  # retrieve level codes
  z = categorical(vals)
  c = levelcode.(z)

  # number of values and levels
  n = length(z)
  k = length(levels(z))

  # minimum lag
  δh = minlag

  # initialize matrix with counts
  C = zeros(Int, k, k)

  # update counts along trajectory
  @inbounds for i in 1:(n - 1)
    c₁ = c[i]
    c₂ = c[i + 1]
    h = evaluate(Euclidean(), centroid(dom, i), centroid(dom, i + 1))
    m = ceil(Int, h / 2δh)
    C[c₁, c₁] += m
    C[c₂, c₂] += m
    C[c₁, c₂] += 1
  end

  # add head of first element
  c₁ = c[1]
  h = evaluate(Euclidean(), centroid(dom, 1), centroid(dom, 2))
  m = ceil(Int, h / 2δh)
  C[c₁, c₁] += m

  # add tail of last element
  c₂ = c[n]
  h = evaluate(Euclidean(), centroid(dom, n - 1), centroid(dom, n))
  m = ceil(Int, h / 2δh)
  C[c₂, c₂] += m

  C
end

# ------------
# THEORETICAL
# ------------

"""
    baseratematrix(l, p)

Transition rate matrix from mean lengths `l` and proportions `p`
that assumes random transitions based on relative proportions.

The transition rate for the pair `i -> j` is given by `-1 / l[i]`
if `i == j` and by `(p[j] / (1 - p[i])) / l[i]` otherwise.

## References

* Carle et al 1998. [Conditional Simulation of Hydrofacies Architecture:
  A Transition Probability/Markov Approach](https://doi.org/10.2110/sepmcheg.01.147)
"""
function baseratematrix(l, p)
  nₗ = length(l)
  nₚ = length(p)

  # sanity checks
  if nₗ != nₚ
    throw(ArgumentError("lengths and proportions must have the same length"))
  end
  if !all(pᵢ -> 0 ≤ pᵢ ≤ 1, p)
    throw(ArgumentError("proportions must be in interval [0, 1]"))
  end

  # Eq. 17 and 18 of Carle et al 1998.
  map(Iterators.product(1:nₗ, 1:nₗ)) do (i, j)
    if i == j
      -1 / l[i]
    else
      (p[j] / (1 - p[i])) / l[i]
    end
  end
end

# -----------------
# HELPER FUNCTIONS
# -----------------

function _domvals(data, var)
  dom = domain(data)
  tab = values(data)
  cols = Tables.columns(tab)
  vals = Tables.getcolumn(cols, Symbol(var))
  dom, vals
end

function _defaultminlag(data)
  d = domain(data)
  n = nelements(d)
  h(i) = evaluate(Euclidean(), centroid(d, i), centroid(d, i + 1))
  minimum(h(i) for i in 1:(n - 1)) / 2
end
