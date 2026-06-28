# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalTransiogram(geotable, var; [options])

Computes the empirical (a.k.a. experimental) omnidirectional
transiogram for categorical variable `var` stored in the `geotable`.

## Options

  * nlags     - number of lags (default to `20`)
  * maxlag    - maximum lag in length units (default to 1/2 of minimum side of bounding box)
  * distance  - custom distance function (default to `Euclidean` distance)
  * lagsearch - lag search method (default to `:ball`)

Available lag search methods:

  * `:full` - loop over all pairs of points available
  * `:ball` - loop over all points within maximum lag

All implemented lag search methods produce the exact same result.
The `:ball` method is considerably faster when the maximum lag is
much smaller than the bounding box of the domain.

See also: [`DirectionalTransiogram`](@ref), [`PlanarTransiogram`](@ref),
[`EmpiricalTransiogramSurface`](@ref).

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)

* Carle et al 1998. [Conditional Simulation of Hydrofacies Architecture:
  A Transition Probability/Markov Approach](https://doi.org/10.2110/sepmcheg.01.147)
"""
struct EmpiricalTransiogram{ℒ<:Len,V,D,E} <: EmpiricalGeoStatsFunction
  counts::Vector{Int}
  abscissas::Vector{ℒ}
  ordinates::Matrix{Vector{V}}
  headcounts::Matrix{Vector{Int}}
  distance::D
  estimator::E
end

function EmpiricalTransiogram(
  data::AbstractGeoTable,
  var;
  nlags=20,
  maxlag=defaultmaxlag(data),
  distance=Euclidean(),
  lagsearch=:ball
)
  # retrieve table and domain
  tab = values(data)
  dom = domain(data)

  # empirical estimators are defined on point sets
  pset = georef(tab, [centroid(dom, i) for i in 1:nelements(dom)])

  # transiograms are estimated based on indicators
  itb = pset |> Select(var) |> OneHot(var)

  # pairs of indicator variables
  ivars = itb |> values |> Tables.columns |> Tables.columnnames
  pairs = Iterators.product(ivars, ivars) |> collect

  # define transiogram estimator
  estim = CarleEstimator()

  # define lag search method
  lsearch = lagsearchmethod(dom, nlags, maxlag, distance, Symbol(lagsearch))

  # perform estimation
  counts, abscissas, ordinates, headcounts = accumulate(itb, pairs, estim, lsearch)

  EmpiricalTransiogram(counts, abscissas, ordinates, headcounts, distance, estim)
end

"""
    EmpiricalTransiogram(partition, var; [options])

Compute the empirical transiogram of the geospatial `partition` for
the categorical variable `var`.

Forwards `options` to the underlying [`EmpiricalTransiogram`](@ref)
calls with geospatial data.
"""
function EmpiricalTransiogram(partition::Partition, var; kwargs...)
  # categorical levels across subsets
  gtb = parent(partition)
  cols = Tables.columns(values(gtb))
  vals = Tables.getcolumn(cols, Symbol(var))
  levs = levels(vals)
  # retain geospatial data with at least two elements
  filtered = Iterators.filter(d -> nelements(domain(d)) > 1, partition)
  @assert !isempty(filtered) "invalid partition of geospatial data"
  t(d) = EmpiricalTransiogram(d |> Levels(var => levs), var; kwargs...)
  tmapreduce(t, merge, collect(filtered))
end

"""
    DirectionalTransiogram(direction, data, var; dtol=1e-6u"m", [options])

Computes the empirical transiogram for the categorical variable `var` stored in
geospatial `data` along a given `direction` with band tolerance `dtol` in length units.

Forwards `options` to the underlying [`EmpiricalTransiogram`](@ref).
"""
function DirectionalTransiogram(dir, data::AbstractGeoTable, var; dtol=1e-6u"m", kwargs...)
  Π = partition(Xoshiro(123), data, DirectionPartition(dir; tol=dtol))
  EmpiricalTransiogram(Π, var; kwargs...)
end

"""
    PlanarTransiogram(normal, data, var; ntol=1e-6u"m", [options])

Computes the empirical transiogram for the categorical variable `var` stored in
geospatial `data` along a plane perpendicular to a `normal` direction with plane
tolerance `ntol` in length units.

Forwards `options` to the underlying [`EmpiricalTransiogram`](@ref).
"""
function PlanarTransiogram(normal, data::AbstractGeoTable, var; ntol=1e-6u"m", kwargs...)
  Π = partition(Xoshiro(123), data, PlanePartition(normal; tol=ntol))
  EmpiricalTransiogram(Π, var; kwargs...)
end

nvariates(t::EmpiricalTransiogram) = size(t.ordinates, 1)

"""
    merge(tα, tβ)

Merge the empirical transiogram `tα` with the empirical transiogram `tβ`
assuming that both transiograms have the same number of lags, distance
and estimator.
"""
function merge(tα::EmpiricalTransiogram{ℒ,V,D,E}, tβ::EmpiricalTransiogram{ℒ,V,D,E}) where {ℒ,V,D,E}
  nα = tα.counts
  nβ = tβ.counts
  xα = tα.abscissas
  xβ = tβ.abscissas
  Yα = tα.ordinates
  Yβ = tβ.ordinates
  Cα = tα.headcounts
  Cβ = tβ.headcounts

  # copy distance and estimator
  d = tα.distance
  e = tα.estimator

  # merge coordinates and bin counts
  n = nα + nβ
  x = @. (xα * nα + xβ * nβ) / n
  Y = map(zip(Yα, Cα, Yβ, Cβ)) do (yα, cα, yβ, cβ)
    @. (yα * cα + yβ * cβ) / (cα + cβ)
  end
  C = map(zip(Cα, Cβ)) do (cα, cβ)
    @. cα + cβ
  end

  # adjust empty bins
  x[n .== 0] .= xα[n .== 0]
  for (y, c) in zip(Y, C)
    y[c .== 0] .= 0
  end

  EmpiricalTransiogram(n, x, Y, C, d, e)
end
