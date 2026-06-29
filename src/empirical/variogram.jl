# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalVariogram(geotable, vars; [options])

Computes the empirical (a.k.a. experimental) omnidirectional
(cross-)variogram for variables `vars` stored in the `geotable`.

## Options

  * nlags     - number of lags (default to `20`)
  * maxlag    - maximum lag in length units (default to 1/2 of minimum side of bounding box)
  * distance  - custom distance function (default to `Euclidean` distance)
  * estimator - variogram estimator (default to `:matheron` estimator)
  * lagsearch - lag search method (default to `:ball`)

Available estimators:

  * `:matheron` - simple estimator based on squared differences
  * `:cressie`  - robust estimator based on 4th power of differences

Available lag search methods:

  * `:full` - loop over all pairs of points available
  * `:ball` - loop over all points within maximum lag

All implemented lag search methods produce the exact same result.
The `:ball` method is considerably faster when the maximum lag is
much smaller than the bounding box of the domain.

See also: [`DirectionalVariogram`](@ref), [`PlanarVariogram`](@ref),
[`EmpiricalVariogramSurface`](@ref).

## References

* Chilès, JP and Delfiner, P. 2012. [Geostatistics: Modeling Spatial Uncertainty]
  (https://onlinelibrary.wiley.com/doi/book/10.1002/9781118136188)

* Webster, R and Oliver, MA. 2007. [Geostatistics for Environmental Scientists]
  (https://onlinelibrary.wiley.com/doi/book/10.1002/9780470517277)

* Hoffimann, J and Zadrozny, B. 2019. [Efficient variography with partition variograms]
  (https://www.sciencedirect.com/science/article/pii/S0098300419302936)
"""
struct EmpiricalVariogram{ℒ<:Len,V,D,E} <: EmpiricalGeoStatsFunction
  counts::Vector{Int}
  abscissas::Vector{ℒ}
  ordinates::Matrix{Vector{V}}
  distance::D
  estimator::E
  variables::Vector{Symbol}
end

function EmpiricalVariogram(
  data::AbstractGeoTable,
  vars;
  nlags=20,
  maxlag=defaultmaxlag(data),
  distance=Euclidean(),
  estimator=:matheron,
  lagsearch=:ball
)
  # selected variables
  sdata = data |> Select(vars)

  # define variogram estimator
  estim = if Symbol(estimator) == :matheron
    MatheronEstimator()
  elseif Symbol(estimator) == :cressie
    CressieEstimator()
  else
    throw(ArgumentError("invalid estimator"))
  end

  # define lag search method
  lsearch = lagsearchmethod(domain(sdata), nlags, maxlag, distance, Symbol(lagsearch))

  # perform estimation
  counts, abscissas, ordinates = _variogram(sdata, estim, lsearch)

  # extract variable names
  names = sdata |> values |> Tables.columns |> Tables.columnnames |> collect

  EmpiricalVariogram(counts, abscissas, ordinates, distance, estim, names)
end

"""
    EmpiricalVariogram(partition, vars; [options])

Computes the empirical (cross-)variogram of the geospatial `partition`
for variables `vars` as described in Hoffimann & Zadrozny 2019.

## References

* Hoffimann, J and Zadrozny, B. 2019. [Efficient variography with partition variograms]
  (https://www.sciencedirect.com/science/article/pii/S0098300419302936)
"""
function EmpiricalVariogram(partition::Partition, vars; kwargs...)
  # retain geospatial data with at least two elements
  filtered = Iterators.filter(d -> nelements(domain(d)) > 1, partition)
  @assert !isempty(filtered) "invalid partition of geospatial data"
  γ(d) = EmpiricalVariogram(d, vars; kwargs...)
  tmapreduce(γ, merge, collect(filtered))
end

issymmetric(::Type{<:EmpiricalVariogram}) = true

nvariables(g::EmpiricalVariogram) = length(g.variables)

variables(g::EmpiricalVariogram) = g.variables

"""
    merge(γα, γβ)

Merge the empirical variogram `γα` with the empirical variogram `γβ`
assuming that both variograms have the same number of lags, distance
and estimator.
"""
function merge(γα::EmpiricalVariogram{ℒ,V,D,E}, γβ::EmpiricalVariogram{ℒ,V,D,E}) where {ℒ,V,D,E}
  nα = γα.counts
  nβ = γβ.counts
  xα = γα.abscissas
  xβ = γβ.abscissas
  Yα = γα.ordinates
  Yβ = γβ.ordinates
  vα = γα.variables
  vβ = γβ.variables

  # copy distance and estimator
  d = γα.distance
  e = γα.estimator

  # merge function for estimator
  mergefun(yα, nα, yβ, nβ) = mergerule(e, yα, nα, yβ, nβ)

  # merge coordinates and bin counts
  n = nα + nβ
  x = @. (xα * nα + xβ * nβ) / n
  Y = map(Yα, Yβ) do yα, yβ
    mergefun.(yα, nα, yβ, nβ)
  end

  # adjust empty bins
  x[n .== 0] .= xα[n .== 0]
  for y in Y
    y[n .== 0] .= 0
  end

  EmpiricalVariogram(n, x, Y, d, e, vα)
end

Base.getindex(g::EmpiricalVariogram, inds::AbstractVector) =
  EmpiricalVariogram(g.counts, g.abscissas, g.ordinates[inds, inds], g.distance, g.estimator, g.variables[inds])

Base.getindex(g::EmpiricalVariogram, ind::Int) =
  EmpiricalVariogram(g.counts, g.abscissas, g.ordinates[[ind], [ind]], g.distance, g.estimator, g.variables[[ind]])

# -------------------------
# CONVENIENCE CONSTRUCTORS
# -------------------------

"""
    DirectionalVariogram(direction, geotable, vars; dtol=1e-6u"m", [options])

Computes the empirical (cross-)variogram for the variables `vars` stored in
the `geotable` along a given `direction` with band tolerance `dtol` in length
units. Forwards `options` to the underlying [`EmpiricalVariogram`](@ref).
"""
function DirectionalVariogram(dir, data::AbstractGeoTable, vars; dtol=1e-6u"m", kwargs...)
  Π = partition(Xoshiro(123), data, DirectionPartition(dir; tol=dtol))
  EmpiricalVariogram(Π, vars; kwargs...)
end

"""
    PlanarVariogram(normal, geotable, vars; ntol=1e-6u"m", [options])

Computes the empirical (cross-)variogram for the variables `vars` stored in
the `geotable` along a plane perpendicular to a `normal` direction with plane
tolerance `ntol` in length units. Forwards `options` to the underlying
[`EmpiricalVariogram`](@ref).
"""
function PlanarVariogram(normal, data::AbstractGeoTable, vars; ntol=1e-6u"m", kwargs...)
  Π = partition(Xoshiro(123), data, PlanePartition(normal; tol=ntol))
  EmpiricalVariogram(Π, vars; kwargs...)
end

# -----------------
# HELPER FUNCTIONS
# -----------------

function _variogram(geotable, estim::VariogramEstimator, lagsearch::LagSearchMethod)
  # lag search parameters
  nlags, maxlag, distance = params(lagsearch)

  # compute lag size
  δh = maxlag / nlags

  # table and domain
  dom = domain(geotable)
  tab = values(geotable)

  # estimators are defined on point sets
  pset = PointSet([centroid(dom, i) for i in 1:nelements(dom)])

  # table columns
  cols = Tables.columns(tab)

  # pairs of variable names
  vars = Tables.columnnames(cols)
  pairs = [(var₁, var₂) for var₁ in vars, var₂ in vars]

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

  # ordinate sums
  Σy = map(pairs) do (var₁, var₂)
    V = returntype(estim, get(var₁), get(var₂))
    zeros(V, nlags)
  end

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
        ns[lag] += 1
        Σx[lag] += h
        for idx in CartesianIndices(pairs)
          # variogram is symmetric
          idx[1] < idx[2] && continue

          # retrieve variables
          var₁, var₂ = pairs[idx]
          z₁ = get(var₁)
          z₂ = get(var₂)

          # evaluate function estimator
          v = accumterm(estim, z₁[i], z₁[j], z₂[i], z₂[j])

          # accumulate if value is valid
          ismissing(v) || (Σy[idx][lag] += v)
        end
      end
    end
  end

  # copy symmetric values to upper triangle
  for idx in CartesianIndices(pairs)
    idx[1] < idx[2] && (Σy[idx] .= Σy[idx[2], idx[1]])
  end

  # bin (or lag) size
  lags = range(δh / 2, stop=maxlag - δh / 2, length=nlags)

  # abscissa
  xs = @. Σx / ns
  xs[ns .== 0] .= lags[ns .== 0]

  # ordinate function
  ordfun(Σ, n) = accumnorm(estim, Σ, n)

  # ordinate
  Y = map(enumerate(pairs)) do (k, _)
    ys = ordfun.(Σy[k], ns)
    ys[ns .== 0] .= zero(eltype(ys))
    ys
  end

  ns, xs, Y
end
