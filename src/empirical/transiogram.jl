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
  * maxlag    - maximum lag in length units (default to 1/10 of minimum side of bounding box)
  * distance  - custom distance function (default to `Euclidean` distance)
  * algorithm - accumulation algorithm (default to `:ball`)

Available algorithms:

  * `:full` - loop over all pairs of points in the data
  * `:ball` - loop over all points inside maximum lag ball

All implemented algorithms produce the exact same result.
The `:ball` algorithm is considerably faster when the
maximum lag is much smaller than the bounding box of
the domain of the data.

See also: [`DirectionalTransiogram`](@ref), [`PlanarTransiogram`](@ref).

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)

* Carle et al 1998. [Conditional Simulation of Hydrofacies Architecture:
  A Transition Probability/Markov Approach](https://doi.org/10.2110/sepmcheg.01.147)
"""
struct EmpiricalTransiogram{ℒ<:Len,V,D}
  counts::Vector{Int}
  abscissa::Vector{ℒ}
  ordinate::Matrix{Vector{V}}
  distance::D
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
  counts, abscissa, ordinate = accumulate(ℐ, pairs, estim, algo)

  EmpiricalTransiogram(counts, abscissa, ordinate, distance)
end
