# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LagBallSearch(nlags, maxlag, distance)

Search pairs of points in geospatial data with
nearest neighbors inside metric ball.
"""
struct LagBallSearch{ℒ<:Len,D} <: LagSearchMethod
  nlags::Int
  maxlag::ℒ
  distance::D
  LagBallSearch(nlags, maxlag::ℒ, distance::D) where {ℒ<:Len,D} = new{float(ℒ),D}(nlags, maxlag, distance)
end

LagBallSearch(nlags, maxlag, distance) = LagBallSearch(nlags, aslen(maxlag), distance)

function neighfun(method::LagBallSearch, pset)
  ball = MetricBall(method.maxlag, method.distance)
  searcher = BallSearch(pset, ball)
  j -> @inbounds(search(pset[j], searcher))
end

skipfun(::LagBallSearch) = (i, j) -> i ≤ j

exitfun(::LagBallSearch) = h -> false
