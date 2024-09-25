# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ratematrix(data, var; [parameters])

Computes the transition rate matrix for categorical variable
`var` stored in geospatial `data` using a minimum lag `minlag`.

See also [`probmatrix`](@ref), [`countmatrix`](@ref)
"""
ratematrix(data::AbstractGeoTable, var; minlag=_defaultminlag(data)) = ratematrix(_domvals(data, var)..., minlag)

"""
    probmatrix(data, var; [parameters])

Computes the transition probability matrix for categorical variable
`var` stored in geospatial `data` using a minimum lag `minlag`.

See also [`ratematrix`](@ref), [`countmatrix`](@ref)
"""
probmatrix(data::AbstractGeoTable, var; minlag=_defaultminlag(data)) = probmatrix(_domvals(data, var)..., minlag)

"""
    countmatrix(data, var; [parameters])

Computes the transition count matrix for categorical variable
`var` stored in geospatial `data` using a minimum lag `minlag`.

See also [`ratematrix`](@ref), [`probmatrix`](@ref)
"""
countmatrix(data::AbstractGeoTable, var; minlag=_defaultminlag(data)) = countmatrix(_domvals(data, var)..., minlag)

function ratematrix(dom, vals, minlag)
  # transition probabilities
  T = probmatrix(dom, vals, minlag)

  # convert probabilities into rates
  Δh = minlag
  log(T) / Δh
end

function probmatrix(dom, vals, minlag)
  # transition counts
  C = countmatrix(dom, vals, minlag)

  # convert counts into probabilities
  C ./ sum(C, dims=2)
end

function countmatrix(dom, vals, minlag)
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
  Δh = minlag

  # initialize matrix with counts
  C = zeros(Int, k, k)

  # update counts along trajectory
  @inbounds for i in 1:(n - 1)
    c₁ = c[i]
    c₂ = c[i + 1]
    h = evaluate(Euclidean(), centroid(dom, i), centroid(dom, i + 1))
    m = ceil(Int, h / 2Δh)
    C[c₁, c₁] += m
    C[c₂, c₂] += m
    C[c₁, c₂] += 1
  end

  # add head of first element
  c₁ = c[1]
  h = evaluate(Euclidean(), centroid(dom, 1), centroid(dom, 2))
  m = ceil(Int, h / 2Δh)
  C[c₁, c₁] += m

  # add tail of last element
  c₂ = c[n]
  h = evaluate(Euclidean(), centroid(dom, n - 1), centroid(dom, n))
  m = ceil(Int, h / 2Δh)
  C[c₂, c₂] += m

  C
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
