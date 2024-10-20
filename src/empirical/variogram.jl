# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalVariogram(data, var₁, var₂=var₁; [parameters])

Computes the empirical (a.k.a. experimental) omnidirectional
(cross-)variogram for variables `var₁` and `var₂` stored in
geospatial `data`.

## Parameters

  * nlags     - number of lags (default to `20`)
  * maxlag    - maximum lag in length units (default to 1/2 of minimum side of bounding box)
  * distance  - custom distance function (default to `Euclidean` distance)
  * estimator - variogram estimator (default to `:matheron` estimator)
  * algorithm - accumulation algorithm (default to `:ball`)

Available estimators:

  * `:matheron` - simple estimator based on squared differences
  * `:cressie`  - robust estimator based on 4th power of differences

Available algorithms:

  * `:full` - loop over all pairs of points in the data
  * `:ball` - loop over all points inside maximum lag ball

All implemented algorithms produce the exact same result.
The `:ball` algorithm is considerably faster when the
maximum lag is much smaller than the bounding box of
the domain of the data.

See also: [`DirectionalVariogram`](@ref), [`PlanarVariogram`](@ref),
[`EmpiricalVarioplane`](@ref).

## References

* Chilès, JP and Delfiner, P. 2012. [Geostatistics: Modeling Spatial Uncertainty]
  (https://onlinelibrary.wiley.com/doi/book/10.1002/9781118136188)

* Webster, R and Oliver, MA. 2007. [Geostatistics for Environmental Scientists]
  (https://onlinelibrary.wiley.com/doi/book/10.1002/9780470517277)

* Hoffimann, J and Zadrozny, B. 2019. [Efficient variography with partition variograms]
  (https://www.sciencedirect.com/science/article/pii/S0098300419302936)
"""
struct EmpiricalVariogram{ℒ<:Len,V,D,E} <: EmpiricalFunction
  counts::Vector{Int}
  abscissas::Vector{ℒ}
  ordinates::Vector{V}
  distance::D
  estimator::E
end

function EmpiricalVariogram(
  data::AbstractGeoTable,
  var₁,
  var₂=var₁;
  nlags=20,
  maxlag=defaultmaxlag(data),
  distance=Euclidean(),
  estimator=:matheron,
  algorithm=:ball
)
  # retrieve table and domain
  𝒯 = values(data)
  𝒟 = domain(data)

  # empirical estimators are defined on point sets
  𝒮 = georef(𝒯, [centroid(𝒟, i) for i in 1:nelements(𝒟)])

  # retrieve estimator and algorithm
  estim, algo = estimalgo(𝒟, nlags, maxlag, distance, estimator, algorithm)

  # accumulate data with chosen algorithm
  counts, abscissas, ordinates = accumulate(𝒮, (var₁, var₂), estim, algo)

  EmpiricalVariogram(counts, abscissas, ordinates, distance, estim)
end

"""
    EmpiricalVariogram(partition, var₁, var₂=var₁; [parameters])

Compute the empirical (cross-)variogram of the geospatial `partition` for
variables `var₁` and `var₂` as described in Hoffimann & Zadrozny 2019.

Optionally, forward `parameters` for the underlying [`EmpiricalVariogram`](@ref).

## References

* Hoffimann, J and Zadrozny, B. 2019. [Efficient variography with partition variograms]
  (https://www.sciencedirect.com/science/article/pii/S0098300419302936)
"""
function EmpiricalVariogram(partition::Partition, var₁, var₂=var₁; kwargs...)
  # retain geospatial data with at least two elements
  filtered = Iterators.filter(d -> nelements(domain(d)) > 1, partition)
  @assert !isempty(filtered) "invalid partition of geospatial data"
  γ(d) = EmpiricalVariogram(d, var₁, var₂; kwargs...)
  tmapreduce(γ, merge, collect(filtered))
end

"""
    DirectionalVariogram(direction, data, var₁, var₂=var₁; dtol=1e-6u"m", [parameters])

Computes the empirical (cross-)variogram for the variables `var₁` and `var₂` stored in
geospatial `data` along a given `direction` with band tolerance `dtol` in length units.

Optionally, forward `parameters` for the underlying [`EmpiricalVariogram`](@ref).
"""
function DirectionalVariogram(dir, data::AbstractGeoTable, var₁, var₂=var₁; dtol=1e-6u"m", kwargs...)
  Π = partition(MersenneTwister(123), data, DirectionPartition(dir; tol=dtol))
  EmpiricalVariogram(Π, var₁, var₂; kwargs...)
end

"""
    PlanarVariogram(normal, data, var₁, var₂=var₁; ntol=1e-6u"m", [parameters])

Computes the empirical (cross-)variogram for the variables `var₁` and `var₂` stored in
geospatial `data` along a plane perpendicular to a `normal` direction with plane
tolerance `ntol` in length units.

Optionally, forward `parameters` for the underlying [`EmpiricalVariogram`](@ref).
"""
function PlanarVariogram(normal, data::AbstractGeoTable, var₁, var₂=var₁; ntol=1e-6u"m", kwargs...)
  Π = partition(MersenneTwister(123), data, PlanePartition(normal; tol=ntol))
  EmpiricalVariogram(Π, var₁, var₂; kwargs...)
end

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
  yα = γα.ordinates
  yβ = γβ.ordinates

  # copy distance and estimator
  d = γα.distance
  e = γα.estimator

  # merge function for estimator
  mergefun(yα, nα, yβ, nβ) = combine(e, yα, nα, yβ, nβ)

  # merge coordinates and bin counts
  n = nα + nβ
  x = @. (xα * nα + xβ * nβ) / n
  y = @. mergefun(yα, nα, yβ, nβ)

  # adjust empty bins
  x[n .== 0] .= xα[n .== 0]
  y[n .== 0] .= 0

  EmpiricalVariogram(n, x, y, d, e)
end
