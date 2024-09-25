# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    AccumAlgorithm

Algorithm for accumulating pairs of points in the
estimation of geostatistical functions.
"""
abstract type AccumAlgorithm end

"""
    accumulate(data, var₁, var₂, estimator, algo)

Accumulate pairs of points in `data` for variables
`var₁` and `var₂` with `estimator` and accumulation
algorithm `algo`.
"""
function accumulate(data, var₁, var₂, estimator::Estimator, algo::AccumAlgorithm)
  # retrieve algorithm parameters
  nlags = algo.nlags
  maxlag = algo.maxlag
  distance = algo.distance

  # compute lag size
  δh = maxlag / nlags

  # table and point set
  𝒯 = values(data)
  𝒫 = domain(data)

  # vectors for variables
  cols = Tables.columns(𝒯)
  z₁ = Tables.getcolumn(cols, Symbol(var₁))
  z₂ = Tables.getcolumn(cols, Symbol(var₂))

  # neighbors function
  neighbors = neighfun(algo, 𝒫)

  # skip condition
  skip = skipfun(algo)

  # early exit condition
  exit = exitfun(algo)

  # accumulation type
  V = returntype(estimator, z₁, z₂)

  # lag sums and counts
  ℒ = Meshes.lentype(𝒫)
  Σx = zeros(ℒ, nlags)
  Σy = zeros(V, nlags)
  ns = zeros(Int, nlags)

  # loop over points inside ball
  @inbounds for j in 1:nelements(𝒫)
    pⱼ = 𝒫[j]
    z₁ⱼ = z₁[j]
    z₂ⱼ = z₂[j]
    for i in neighbors(j)
      # skip to avoid double counting
      skip(i, j) && continue

      pᵢ = 𝒫[i]
      z₁ᵢ = z₁[i]
      z₂ᵢ = z₂[i]

      # evaluate geospatial lag
      h = evaluate(distance, pᵢ, pⱼ)

      # early exit if out of range
      exit(h) && continue

      # evaluate (cross-)variance
      v = formula(estimator, z₁ᵢ, z₁ⱼ, z₂ᵢ, z₂ⱼ)

      # bin (or lag) where to accumulate result
      lag = ceil(Int, h / δh)
      lag == 0 && @warn "duplicate coordinates found, consider using `UniqueCoords`"

      if 0 < lag ≤ nlags && !ismissing(v)
        Σx[lag] += h
        Σy[lag] += v
        ns[lag] += 1
      end
    end
  end

  # bin (or lag) size
  lags = range(δh / 2, stop=maxlag - δh / 2, length=nlags)

  # ordinate function
  ordfun(Σy, n) = normsum(estimator, Σy, n)

  # variogram abscissa
  xs = @. Σx / ns
  xs[ns .== 0] .= lags[ns .== 0]

  # variogram ordinate
  ys = @. ordfun(Σy, ns)
  ys[ns .== 0] .= zero(eltype(ys))

  xs, ys, ns
end

include("algorithms/fullsearch.jl")
include("algorithms/ballsearch.jl")
