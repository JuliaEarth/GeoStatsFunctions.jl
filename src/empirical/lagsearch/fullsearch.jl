# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LagFullSearch(nlags, maxlag, distance)

Search pairs of points in geospatial data with
exhaustive (or full) search.
"""
struct LagFullSearch{ℒ<:Len,D} <: LagSearchMethod
  nlags::Int
  maxlag::ℒ
  distance::D
  LagFullSearch(nlags, maxlag::ℒ, distance::D) where {ℒ<:Len,D} = new{float(ℒ),D}(nlags, maxlag, distance)
end

LagFullSearch(nlags, maxlag, distance) = LagFullSearch(nlags, aslen(maxlag), distance)

neighfun(::LagFullSearch, pset) = j -> (j + 1):nelements(pset)

skipfun(::LagFullSearch) = (i, j) -> false

exitfun(method::LagFullSearch) = h -> h > method.maxlag
