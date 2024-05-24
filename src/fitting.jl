# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# models that can be fitted currently
fittable() = filter(isstationary, setdiff(subtypes(Variogram), (NuggetEffect, NestedVariogram)))

"""
    VariogramFitAlgo

An algorithm for fitting theoretical variograms.
"""
abstract type VariogramFitAlgo end

"""
    WeightedLeastSquares()
    WeightedLeastSquares(w)

Fit theoretical variogram using weighted least squares with weighting
function `w` (e.g. h -> 1/h). If no weighting function is provided,
bin counts of empirical variogram are normalized and used as weights.
"""
struct WeightedLeastSquares <: VariogramFitAlgo
  weightfun::Union{Function,Nothing}
end

WeightedLeastSquares() = WeightedLeastSquares(nothing)

"""
    fit(V, g, algo=WeightedLeastSquares(); range=nothing, sill=nothing, nugget=nothing)

Fit theoretical variogram type `V` to empirical variogram `g`
using algorithm `algo`.

Optionally fix `range`, `sill` or `nugget` by passing them as keyword arguments, or
set their maximum value with `maxrange`, `maxsill` or `maxnugget`.

## Examples

```julia
julia> fit(SphericalVariogram, g)
julia> fit(ExponentialVariogram, g)
julia> fit(ExponentialVariogram, g, sill=1.0)
julia> fit(ExponentialVariogram, g, maxsill=1.0)
julia> fit(GaussianVariogram, g, WeightedLeastSquares())
```
"""
fit(V::Type{<:Variogram}, g::EmpiricalVariogram, algo::VariogramFitAlgo=WeightedLeastSquares(); kwargs...) =
  fit_impl(V, g, algo; kwargs...) |> first

"""
    fit(Vs, g, algo=WeightedLeastSquares(); kwargs...)

Fit theoretical variogram types `Vs` to empirical variogram `g`
using algorithm `algo` and return the one with minimum error.

## Examples

```julia
julia> fit([SphericalVariogram, ExponentialVariogram], g)
```
"""
function fit(Vs, g::EmpiricalVariogram, algo::VariogramFitAlgo=WeightedLeastSquares(); kwargs...)
  # fit each variogram type
  res = [fit_impl(V, g, algo; kwargs...) for V in Vs]
  γs, ϵs = first.(res), last.(res)

  # return best candidate
  γs[argmin(ϵs)]
end

"""
    fit(Variogram, g, algo=WeightedLeastSquares(); kwargs...)

Fit all "fittable" subtypes of `Variogram` to empirical variogram `g`
using algorithm `algo` and return the one with minimum error.

## Examples

```julia
julia> fit(Variogram, g)
julia> fit(Variogram, g, WeightedLeastSquares())
```

See also `GeoStatsFunctions.fittable()`.
"""
fit(::Type{Variogram}, g::EmpiricalVariogram, algo::VariogramFitAlgo=WeightedLeastSquares(); kwargs...) =
  fit(fittable(), g, algo; kwargs...)

"""
    fit(V, g, weightfun; kwargs...)
    fit(Variogram, g, weightfun; kwargs...)

Convenience method that forwards the weighting function `weightfun`
to the `WeightedLeastSquares` algorithm.

## Examples

```julia
fit(SphericalVariogram, g, h -> exp(-h))
fit(Variogram, g, h -> exp(-h/100))
```
"""
fit(V, g::EmpiricalVariogram, weightfun::Function; kwargs...) = fit(V, g, WeightedLeastSquares(weightfun); kwargs...)

# ---------------
# IMPLEMENTATION
# ---------------

function fit_impl(
  V::Type{<:Variogram},
  g::EmpiricalVariogram,
  algo::WeightedLeastSquares;
  range=nothing,
  sill=nothing,
  nugget=nothing,
  maxrange=nothing,
  maxsill=nothing,
  maxnugget=nothing
)
  # values of empirical variogram
  x, y, n = values(g)

  # custom ball of given radius
  ball(r) = MetricBall(r, distance(g))

  # discard invalid bins
  x = x[n .> 0]
  y = y[n .> 0]
  n = n[n .> 0]

  # strip units of coordinates
  ux = unit(eltype(x))
  uy = unit(eltype(y))
  x′ = ustrip.(x)
  y′ = ustrip.(y)

  # strip units of kwargs
  range′ = isnothing(range) ? range : ustrip(ux, range)
  sill′ = isnothing(sill) ? sill : ustrip(uy, sill)
  nugget′ = isnothing(nugget) ? nugget : ustrip(uy, nugget)
  maxrange′ = isnothing(maxrange) ? maxrange : ustrip(ux, maxrange)
  maxsill′ = isnothing(maxsill) ? maxsill : ustrip(uy, maxsill)
  maxnugget′ = isnothing(maxnugget) ? maxnugget : ustrip(uy, maxnugget)

  # evaluate weights
  f = algo.weightfun
  w = isnothing(f) ? n / sum(n) : map(xᵢ -> ustrip(f(xᵢ)), x)

  # objective function
  function J(θ)
    γ = V(ball(θ[1]), sill=θ[2], nugget=θ[3])
    sum(i -> w[i] * (γ(x′[i]) - y′[i])^2, eachindex(w, x′, y′))
  end

  # linear constraint (sill ≥ nugget)
  L(θ) = θ[2] ≥ θ[3] ? 0.0 : θ[3] - θ[2]

  # penalty for linear constraint (J + λL)
  λ = sum(yᵢ -> yᵢ^2, y′)

  # maximum range, sill and nugget
  xmax = maximum(x′)
  ymax = maximum(y′)
  rmax = isnothing(maxrange′) ? xmax : maxrange′
  smax = isnothing(maxsill′) ? ymax : maxsill′
  nmax = isnothing(maxnugget′) ? ymax : maxnugget′

  # initial guess
  rₒ = isnothing(range′) ? rmax / 3 : range′
  sₒ = isnothing(sill′) ? 0.95 * smax : sill′
  nₒ = isnothing(nugget′) ? 1e-6 : nugget′
  θₒ = [rₒ, sₒ, nₒ]

  # box constraints
  δ = 1e-8
  rₗ, rᵤ = isnothing(range′) ? (zero(rmax), rmax) : (range′ - δ, range′ + δ)
  sₗ, sᵤ = isnothing(sill′) ? (zero(smax), smax) : (sill′ - δ, sill′ + δ)
  nₗ, nᵤ = isnothing(nugget′) ? (zero(nmax), nmax) : (nugget′ - δ, nugget′ + δ)
  l = [rₗ, sₗ, nₗ]
  u = [rᵤ, sᵤ, nᵤ]

  # solve optimization problem
  sol = Optim.optimize(θ -> J(θ) + λ * L(θ), l, u, θₒ)
  ϵ = Optim.minimum(sol)
  θ = Optim.minimizer(sol)

  # optimal variogram (with units)
  γ = V(ball(θ[1]), sill=θ[2] * uy, nugget=θ[3] * uy)

  γ, ϵ
end
