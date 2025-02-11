# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalTransiogram(data, var; [parameters])

Computes the empirical (a.k.a. experimental) omnidirectional
transiogram for categorical variable `var` stored in geospatial
`data`.

## Parameters

  * nlags     - number of lags (default to `20`)
  * maxlag    - maximum lag in length units (default to 1/2 of minimum side of bounding box)
  * distance  - custom distance function (default to `Euclidean` distance)
  * algorithm - accumulation algorithm (default to `:ball`)

Available algorithms:

  * `:full` - loop over all pairs of points in the data
  * `:ball` - loop over all points inside maximum lag ball

All implemented algorithms produce the exact same result.
The `:ball` algorithm is considerably faster when the
maximum lag is much smaller than the bounding box of
the domain of the data.

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
  distance::D
  estimator::E
end

function EmpiricalTransiogram(
  data::AbstractGeoTable,
  var;
  nlags=20,
  maxlag=defaultmaxlag(data),
  distance=Euclidean(),
  algorithm=:ball
)
  # retrieve table and domain
  𝒯 = values(data)
  𝒟 = domain(data)

  # empirical estimators are defined on point sets
  𝒮 = georef(𝒯, [centroid(𝒟, i) for i in 1:nelements(𝒟)])

  # transiograms are estimated based on indicators
  ℐ = 𝒮 |> OneHot(var)

  # pairs of indicator variables
  ivars = ℐ |> values |> Tables.columns |> Tables.columnnames
  pairs = Iterators.product(ivars, ivars) |> collect

  # retrieve estimator and algorithm
  estim, algo = estimalgo(𝒟, nlags, maxlag, distance, :carle, algorithm)

  # accumulate data with chosen algorithm
  counts, abscissas, ordinates = accumulate(ℐ, pairs, estim, algo)

  EmpiricalTransiogram(counts, abscissas, ordinates, distance, estim)
end

"""
    EmpiricalTransiogram(partition, var; [parameters])

Compute the empirical transiogram of the geospatial `partition` for
the categorical variable `var`.

Optionally, forward `parameters` for the underlying [`EmpiricalTransiogram`](@ref).
"""
function EmpiricalTransiogram(partition::Partition, var; kwargs...)
  # retain categorical levels across subsets
  levs = levels(parent(partition)[:, var])
  # retain geospatial data with at least two elements
  filtered = Iterators.filter(d -> nelements(domain(d)) > 1, partition)
  @assert !isempty(filtered) "invalid partition of geospatial data"
  t(d) = EmpiricalTransiogram(d |> Levels(var => levs), var; kwargs...)
  tmapreduce(t, merge, collect(filtered))
end

"""
    DirectionalTransiogram(direction, data, var; dtol=1e-6u"m", [parameters])

Computes the empirical transiogram for the categorical variable `var` stored in
geospatial `data` along a given `direction` with band tolerance `dtol` in length units.

Optionally, forward `parameters` for the underlying [`EmpiricalTransiogram`](@ref).
"""
function DirectionalTransiogram(dir, data::AbstractGeoTable, var; dtol=1e-6u"m", kwargs...)
  Π = partition(MersenneTwister(123), data, DirectionPartition(dir; tol=dtol))
  EmpiricalTransiogram(Π, var; kwargs...)
end

"""
    PlanarTransiogram(normal, data, var; ntol=1e-6u"m", [parameters])

Computes the empirical transiogram for the categorical variable `var` stored in
geospatial `data` along a plane perpendicular to a `normal` direction with plane
tolerance `ntol` in length units.

Optionally, forward `parameters` for the underlying [`EmpiricalTransiogram`](@ref).
"""
function PlanarTransiogram(normal, data::AbstractGeoTable, var; ntol=1e-6u"m", kwargs...)
  Π = partition(MersenneTwister(123), data, PlanePartition(normal; tol=ntol))
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

  # copy distance and estimator
  d = tα.distance
  e = tα.estimator

  # merge function for estimator
  mergefun(yα, nα, yβ, nβ) = combine(e, yα, nα, yβ, nβ)

  # merge coordinates and bin counts
  n = nα + nβ
  x = @. (xα * nα + xβ * nβ) / n
  Y = map(zip(Yα, Yβ)) do (yα, yβ)
    @. mergefun(yα, nα, yβ, nβ)
  end

  # adjust empty bins
  x[n .== 0] .= xα[n .== 0]
  for y in Y
    y[n .== 0] .= 0
  end

  EmpiricalTransiogram(n, x, Y, d, e)
end
