# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    EmpiricalTransiogram(counts, abscissas, ordinates, headcounts, distance, estimator, variables)

Type that stores all information about the empirical (a.k.a. experimental) transiogram.
It is returned by the [`transiogram`](@ref) function and is considered an internal type
(i.e., not intended for end-users).
"""
struct EmpiricalTransiogram{ℒ<:Len,V,D,E} <: EmpiricalGeoStatsFunction
  counts::Vector{Int}
  abscissas::Vector{ℒ}
  ordinates::Matrix{Vector{V}}
  headcounts::Matrix{Vector{Int}}
  distance::D
  estimator::E
  variables::Vector{Symbol}
end

issymmetric(::Type{<:EmpiricalTransiogram}) = false

nvariables(t::EmpiricalTransiogram) = length(t.variables)

variables(t::EmpiricalTransiogram) = t.variables

function merge(tα::EmpiricalTransiogram{ℒ,V,D,E}, tβ::EmpiricalTransiogram{ℒ,V,D,E}) where {ℒ,V,D,E}
  nα = tα.counts
  nβ = tβ.counts
  xα = tα.abscissas
  xβ = tβ.abscissas
  Yα = tα.ordinates
  Yβ = tβ.ordinates
  Cα = tα.headcounts
  Cβ = tβ.headcounts
  vα = tα.variables
  vβ = tβ.variables

  # copy distance and estimator
  d = tα.distance
  e = tα.estimator

  # merge coordinates and bin counts
  n = nα + nβ
  x = @. (xα * nα + xβ * nβ) / n
  Y = map(Yα, Cα, Yβ, Cβ) do yα, cα, yβ, cβ
    @. (yα * cα + yβ * cβ) / (cα + cβ)
  end
  C = map(Cα, Cβ) do cα, cβ
    @. cα + cβ
  end

  # adjust empty bins
  x[n .== 0] .= xα[n .== 0]
  for (y, c) in zip(Y, C)
    y[c .== 0] .= 0
  end

  EmpiricalTransiogram(n, x, Y, C, d, e, vα)
end

Base.getindex(t::EmpiricalTransiogram, inds::AbstractVector) = EmpiricalTransiogram(
  t.counts,
  t.abscissas,
  t.ordinates[inds, inds],
  t.headcounts[inds, inds],
  t.distance,
  t.estimator,
  t.variables[inds]
)

Base.getindex(t::EmpiricalTransiogram, ind::Int) = EmpiricalTransiogram(
  t.counts,
  t.abscissas,
  t.ordinates[[ind], [ind]],
  t.headcounts[[ind], [ind]],
  t.distance,
  t.estimator,
  t.variables[[ind]]
)

# -------------------
# END-USER INTERFACE
# -------------------

"""
    transiogram(geotable, [var]; [options])

Computes the empirical (a.k.a. experimental) transiogram for
categorical variable `var` stored in the `geotable`. The variable
can be specified by name or index, and the first variable is used
by default.

    transiogram(partition, [var]; [options])

Alternatively, computes the empirical transiogram of the geospatial
`partition` as described in Hoffimann & Zadrozny 2019.

## Options

  * dir       - direction for directional transiogram (default to `nothing`)
  * dirtol    - tolerance for directional transiogram (default to `0.5u"m"`)
  * maxlag    - maximum lag in length units (default to 1/2 of minimum side of bounding box)
  * nlags     - number of lags (default to `20`)
  * distance  - custom distance function (default to `Euclidean` distance)
  * lagsearch - lag search method (default to `:ball`)

Available lag search methods:

  * `:full` - loop over all pairs of points available
  * `:ball` - loop over all points within maximum lag

All implemented lag search methods produce the exact same result.
The `:ball` method is considerably faster when the maximum lag is
much smaller than the bounding box of the domain.

See also [`transiogramsurface`](@ref).

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)

* Carle et al 1998. [Conditional Simulation of Hydrofacies Architecture:
  A Transition Probability/Markov Approach](https://doi.org/10.2110/sepmcheg.01.147)

* Hoffimann, J and Zadrozny, B. 2019. [Efficient variography with
  partition variograms](https://www.sciencedirect.com/science/article/pii/S0098300419302936)
"""
function transiogram(
  geotable::AbstractGeoTable,
  var=1;
  dir=nothing,
  dirtol=0.5u"m",
  maxlag=defaultmaxlag(geotable),
  nlags=20,
  distance=Euclidean(),
  lagsearch=:ball
)
  if isnothing(dir)
    _transiogram(geotable, var; maxlag, nlags, distance, lagsearch)
  else
    part = partition(Xoshiro(123), geotable, DirectionPartition(dir; tol=dirtol))
    transiogram(part, var; maxlag, nlags, distance, lagsearch)
  end
end

function transiogram(part::Partition, var=1; kwargs...)
  # categorical levels across subsets
  levs = let
    gtb = parent(part) |> Select(var)
    cols = Tables.columns(values(gtb))
    vals = Tables.getcolumn(cols, 1)
    levels(vals)
  end
  # retain geospatial data with at least two elements
  filtered = Iterators.filter(gtb -> nelements(domain(gtb)) > 1, part)
  @assert !isempty(filtered) "invalid partition of geospatial data, try increasing tolerance parameters"
  trans(gtb) = _transiogram(gtb |> Levels(var => levs), var; kwargs...)
  tmapreduce(trans, merge, collect(filtered))
end

# -----------------
# HELPER FUNCTIONS
# -----------------

function _transiogram(geotable, var; maxlag, nlags, distance, lagsearch)
  # indicators of categorical variable
  gtb = geotable |> Select(var) |> OneHot(var)

  # define transiogram estimator
  estim = CarleEstimator()

  # define lag search method
  lsearch = lagsearchmethod(domain(gtb), nlags, maxlag, distance, Symbol(lagsearch))

  # perform estimation
  counts, abscissas, ordinates, headcounts = _transiogramestimate(gtb, estim, lsearch)

  # extract variable names
  names = Symbol.(levels(geotable[:, var]))

  EmpiricalTransiogram(counts, abscissas, ordinates, headcounts, distance, estim, names)
end

function _transiogramestimate(geotable, ::CarleEstimator, lagsearch::LagSearchMethod)
  # lag search parameters
  nlags, maxlag, distance = params(lagsearch)

  # compute lag size
  δh = maxlag / nlags

  # table and domain
  dom = domain(geotable)
  tab = values(geotable)

  # estimators are defined on point sets
  pset = PointSet([centroid(dom, i) for i in 1:nelements(dom)])

  # columns with indicator variables
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
        ns[lag] += 1
        Σx[lag] += h
        for idx in CartesianIndices(pairs)
          # retrieve indicator variables
          var₁, var₂ = pairs[idx]
          I₁ = get(var₁)
          I₂ = get(var₂)

          # accumulate numerator and denominator
          Σn[idx][lag] += I₁[i] * I₂[j]
          Σd[idx][lag] += I₁[i]
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
  Y = map(CartesianIndices(pairs)) do idx
    ys = Σn[idx] ./ Σd[idx]
    ys[Σd[idx] .== 0] .= zero(eltype(ys))
    ys
  end

  # head count
  C = Σd

  ns, xs, Y, C
end
