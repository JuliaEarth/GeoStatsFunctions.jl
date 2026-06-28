# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LagSearchMethod

Algorithm used for accumulating values in
the estimation of geostatistical functions.
"""
abstract type LagSearchMethod end

include("lagsearch/fullsearch.jl")
include("lagsearch/ballsearch.jl")

"""
    lagsearchmethod(domain, nlags, maxlag, distance, option)

Determine valid lag search method from user inputs and given
`option` which can be either `:full` or `:ball`.
"""
function lagsearchmethod(dom, nlags, maxlag, distance, option)
  # sanity checks
  @assert nelements(dom) > 1 "domain must contain at least two points"
  @assert nlags > 0 "number of lags must be positive"
  @assert maxlag > zero(maxlag) "maximum lag must be positive"
  @assert option ∈ (:full, :ball) "invalid lag search method"

  # ball search with NearestNeighbors.jl requires AbstractFloat and MinkowskiMetric
  # https://github.com/KristofferC/NearestNeighbors.jl/issues/13
  isfloat = Unitful.numtype(Meshes.lentype(dom)) <: AbstractFloat
  isminkowski = distance isa MinkowskiMetric

  # warn users requesting :ball option with invalid parameters
  (option == :ball && !isfloat) &&
    @warn "lag ball search requires floating point coordinates, falling back to lag full search"
  (option == :ball && !isminkowski) &&
    @warn "lag ball search requires Minkowski metric, falling back to lag full search"

  # choose lag search method
  if option == :ball && isfloat && isminkowski
    LagBallSearch(nlags, maxlag, distance)
  else
    LagFullSearch(nlags, maxlag, distance)
  end
end
