# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BallSearchAccum(nlags, maxlag, distance)

Accumulate pairs of points in geospatial data with
nearest neighbors inside metric ball.
"""
struct BallSearchAccum{ℒ<:Len,D} <: VariogramAccumAlgo
  nlags::Int
  maxlag::ℒ
  distance::D
  BallSearchAccum(nlags, maxlag::ℒ, distance::D) where {ℒ<:Len,D} = new{float(ℒ),D}(nlags, maxlag, distance)
end

BallSearchAccum(nlags, maxlag, distance) = BallSearchAccum(nlags, addunit(maxlag, u"m"), distance)

function neighfun(algo::BallSearchAccum, pset)
  ball = MetricBall(algo.maxlag, algo.distance)
  searcher = BallSearch(pset, ball)
  j -> @inbounds(search(pset[j], searcher))
end

skipfun(::BallSearchAccum) = (i, j) -> i ≤ j

exitfun(::BallSearchAccum) = h -> false
