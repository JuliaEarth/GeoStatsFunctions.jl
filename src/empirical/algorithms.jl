# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    AccumAlgo

Algorithm used for accumulating values in
the estimation of geostatistical functions.
"""
abstract type AccumAlgo end

"""
    accumulate(data, pairs, estimator, algo)

Accumulate values for `pairs` of variables stored
in `data` with `estimator` and accumulation `algo`.
"""
function accumulate(data, (var₁, var₂), estimator::Estimator, algo::AccumAlgo)
  # retrieve algorithm parameters
  nlags = algo.nlags
  maxlag = algo.maxlag
  distance = algo.distance

  # compute lag size
  δh = maxlag / nlags

  # table and point set
  𝒯 = values(data)
  𝒫 = domain(data)

  # table columns
  cols = Tables.columns(𝒯)

  # get column from variable name
  get(var) = Tables.getcolumn(cols, Symbol(var))

  # vectors for variables
  z₁ = get(var₁)
  z₂ = get(var₂)

  # neighbors function
  neighbors = neighfun(algo, 𝒫)

  # skip condition
  skip = skipfun(algo)

  # early exit condition
  exit = exitfun(algo)

  # lag counts and abscissa sums
  ℒ = Meshes.lentype(𝒫)
  ns = zeros(Int, nlags)
  Σx = zeros(ℒ, nlags)

  # ordinate sums
  V = returntype(estimator, z₁, z₂)
  Σy = zeros(V, nlags)

  # loop over pairs of points
  @inbounds for j in 1:nelements(𝒫)
    pⱼ = 𝒫[j]
    for i in neighbors(j)
      # skip to avoid double counting
      skip(i, j) && continue

      pᵢ = 𝒫[i]

      # evaluate geospatial lag
      h = evaluate(distance, pᵢ, pⱼ)

      # early exit if out of range
      exit(h) && continue

      # bin (or lag) where to accumulate result
      lag = ceil(Int, h / δh)
      lag == 0 && @warn "duplicate coordinates found, consider using `UniqueCoords`"

      # accumulate if lag is valid
      if 0 < lag ≤ nlags
        # evaluate function estimator
        v = accumterm(estimator, z₁[i], z₁[j], z₂[i], z₂[j])

        # accumulate if value is valid
        if !ismissing(v)
          ns[lag] += 1
          Σx[lag] += h
          Σy[lag] += v
        end
      end
    end
  end

  # ordinate function
  ordfun(Σy, n) = accumnorm(estimator, Σy, n)

  # bin (or lag) size
  lags = range(δh / 2, stop=maxlag - δh / 2, length=nlags)

  # abscissa
  xs = @. Σx / ns
  xs[ns .== 0] .= lags[ns .== 0]

  # ordinate
  ys = @. ordfun(Σy, ns)
  ys[ns .== 0] .= zero(eltype(ys))

  ns, xs, ys
end

function accumulate(data, pairs, ::CarleEstimator, algo::AccumAlgo)
  # retrieve algorithm parameters
  nlags = algo.nlags
  maxlag = algo.maxlag
  distance = algo.distance

  # compute lag size
  δh = maxlag / nlags

  # table and point set
  𝒯 = values(data)
  𝒫 = domain(data)

  # table columns
  cols = Tables.columns(𝒯)

  # get column from variable name
  get(var) = Tables.getcolumn(cols, Symbol(var))

  # neighbors function
  neighbors = neighfun(algo, 𝒫)

  # skip condition
  skip = skipfun(algo)

  # early exit condition
  exit = exitfun(algo)

  # lag counts and abscissa sums
  ℒ = Meshes.lentype(𝒫)
  ns = zeros(Int, nlags)
  Σx = zeros(ℒ, nlags)

  # numerator and denominator sums
  Σn = map(_ -> zeros(Int, nlags), pairs)
  Σd = map(_ -> zeros(Int, nlags), pairs)

  # loop over pairs of points
  @inbounds for j in 1:nelements(𝒫)
    pⱼ = 𝒫[j]
    for i in neighbors(j)
      # skip to avoid double counting
      skip(i, j) && continue

      pᵢ = 𝒫[i]

      # evaluate geospatial lag
      h = evaluate(distance, pᵢ, pⱼ)

      # early exit if out of range
      exit(h) && continue

      # bin (or lag) where to accumulate result
      lag = ceil(Int, h / δh)
      lag == 0 && @warn "duplicate coordinates found, consider using `UniqueCoords`"

      # accumulate if lag is valid
      if 0 < lag ≤ nlags
        for (k, (var₁, var₂)) in enumerate(pairs)
          # retrieve indicator variables
          I₁ = get(var₁)
          I₂ = get(var₂)

          # accumulate numerator and denominator
          ns[lag] += 1
          Σx[lag] += h
          Σn[k][lag] += I₁[i] * I₂[j]
          Σd[k][lag] += I₁[i]
        end
      end
    end
  end

  # bin (or lag) size
  lags = range(δh / 2, stop=maxlag - δh / 2, length=nlags)

  # abscissa
  xs = @. Σx / ns
  xs[ns .== 0] .= lags[ns .== 0]

  # ordinate
  Y = map(enumerate(pairs)) do (k, _)
    ys = Σn[k] ./ Σd[k]
    ys[Σd[k] .== 0] .= zero(eltype(ys))
    ys
  end

  # head count
  C = Σd

  ns, xs, Y, C
end

include("algorithms/fullsearch.jl")
include("algorithms/ballsearch.jl")
