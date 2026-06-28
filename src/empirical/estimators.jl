# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FunctionEstimator

An estimator of geostatistical functions.
"""
abstract type FunctionEstimator end

# ---------------------
# VARIOGRAM ESTIMATORS
# ---------------------

"""
    VariogramEstimator

An estimator of variogram functions.
"""
abstract type VariogramEstimator <: FunctionEstimator end

returntype(e::VariogramEstimator, z₁, z₂) = typeof(accumterm(e, z₁[1], z₁[2], z₂[1], z₂[2]))

include("estimators/matheron.jl")
include("estimators/cressie.jl")

function accumulate(geotable, (var₁, var₂), e::VariogramEstimator, lagsearch::LagSearchMethod)
  # lag search parameters
  nlags, maxlag, distance = params(lagsearch)

  # compute lag size
  δh = maxlag / nlags

  # table and domain
  tab = values(geotable)
  dom = domain(geotable)

  # estimators are defined on point sets
  pset = PointSet([centroid(dom, i) for i in 1:nelements(dom)])

  # table columns
  cols = Tables.columns(tab)

  # get column from variable name
  get(var) = Tables.getcolumn(cols, Symbol(var))

  # vectors for variables
  z₁ = get(var₁)
  z₂ = get(var₂)

  # neighbors function
  neighbors = neighfun(lagsearch, pset)

  # skip condition
  skip = skipfun(lagsearch)

  # early exit condition
  exit = exitfun(lagsearch)

  # lag counts and abscissa sums
  ℒ = Meshes.lentype(pset)
  ns = zeros(Int, nlags)
  Σx = zeros(ℒ, nlags)

  # ordinate sums
  V = returntype(e, z₁, z₂)
  Σy = zeros(V, nlags)

  # loop over pairs of points
  @inbounds for j in 1:nelements(pset)
    pⱼ = pset[j]
    for i in neighbors(j)
      # skip to avoid double counting
      skip(i, j) && continue

      pᵢ = pset[i]

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
        v = accumterm(e, z₁[i], z₁[j], z₂[i], z₂[j])

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
  ordfun(Σy, n) = accumnorm(e, Σy, n)

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

# -----------------------
# TRANSIOGRAM ESTIMATORS
# -----------------------

"""
    TransiogramEstimator

An estimator of transiogram functions.
"""
abstract type TransiogramEstimator <: FunctionEstimator end

include("estimators/carle.jl")

function accumulate(geotable, pairs, ::CarleEstimator, lagsearch::LagSearchMethod)
  # lag search parameters
  nlags, maxlag, distance = params(lagsearch)

  # compute lag size
  δh = maxlag / nlags

  # table and domain
  tab = values(geotable)
  dom = domain(geotable)

  # estimators are defined on point sets
  pset = PointSet([centroid(dom, i) for i in 1:nelements(dom)])

  # table columns
  cols = Tables.columns(tab)

  # get column from variable name
  get(var) = Tables.getcolumn(cols, Symbol(var))

  # neighbors function
  neighbors = neighfun(lagsearch, pset)

  # skip condition
  skip = skipfun(lagsearch)

  # early exit condition
  exit = exitfun(lagsearch)

  # lag counts and abscissa sums
  ℒ = Meshes.lentype(pset)
  ns = zeros(Int, nlags)
  Σx = zeros(ℒ, nlags)

  # numerator and denominator sums
  Σn = map(_ -> zeros(Int, nlags), pairs)
  Σd = map(_ -> zeros(Int, nlags), pairs)

  # loop over pairs of points
  @inbounds for j in 1:nelements(pset)
    pⱼ = pset[j]
    for i in neighbors(j)
      # skip to avoid double counting
      skip(i, j) && continue

      pᵢ = pset[i]

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
