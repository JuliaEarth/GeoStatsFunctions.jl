# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    estimalgo(domain, nlags, maxlag, distance, estimator, algorithm)

Retrieve concrete estimator and algorithm from user input.
"""
function estimalgo(dom, nlags, maxlag, distance, estimator, algorithm)
  # sanity checks
  @assert nelements(dom) > 1 "variogram requires at least 2 elements"
  @assert nlags > 0 "number of lags must be positive"
  @assert maxlag > zero(maxlag) "maximum lag must be positive"
  @assert estimator ∈ (:matheron, :cressie) "invalid empirical estimator"
  @assert algorithm ∈ (:full, :ball) "invalid accumulation algorithm"

  # choose empirical estimator
  estim = if estimator == :matheron
    MatheronEstimator()
  elseif estimator == :cressie
    CressieEstimator()
  else
    throw(ArgumentError("invalid estimator"))
  end

  # ball search with NearestNeighbors.jl requires AbstractFloat and MinkowskiMetric
  # https://github.com/KristofferC/NearestNeighbors.jl/issues/13
  isfloat = Unitful.numtype(Meshes.lentype(dom)) <: AbstractFloat
  isminkowski = distance isa MinkowskiMetric

  # warn users requesting :ball option with invalid parameters
  (algorithm == :ball && !isfloat) && @warn ":ball algorithm requires floating point coordinates, falling back to :full"
  (algorithm == :ball && !isminkowski) && @warn ":ball algorithm requires Minkowski metric, falling back to :full"

  # choose accumulation algorithm
  algo = if algorithm == :ball && isfloat && isminkowski
    BallSearchAccum(nlags, maxlag, distance)
  else
    FullSearchAccum(nlags, maxlag, distance)
  end

  estim, algo
end
