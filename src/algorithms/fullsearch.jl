# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FullSearchAccum(nlags, maxlag, distance)

Accumulate pairs of points in geospatial data with
exhaustive (or full) search.
"""
struct FullSearchAccum{ℒ<:Len,D} <: VariogramAccumAlgo
  nlags::Int
  maxlag::ℒ
  distance::D
  FullSearchAccum(nlags, maxlag::ℒ, distance::D) where {ℒ<:Len,D} = new{float(ℒ),D}(nlags, maxlag, distance)
end

FullSearchAccum(nlags, maxlag, distance) = FullSearchAccum(nlags, addunit(maxlag, u"m"), distance)

neighfun(::FullSearchAccum, pset) = j -> (j + 1):nelements(pset)

skipfun(::FullSearchAccum) = (i, j) -> false

exitfun(algo::FullSearchAccum) = h -> h > algo.maxlag
