# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# code to generate list of stationary variogram models:
# filter(isstationary, setdiff(subtypes(Variogram), (NuggetEffect, NestedVariogram))) |> Tuple
stationaryvariograms() = (
  CircularVariogram,
  CubicVariogram,
  ExponentialVariogram,
  GaussianVariogram,
  MaternVariogram,
  PentasphericalVariogram,
  SineHoleVariogram,
  SphericalVariogram
)

"""
    fit(G, g, algo=WeightedLeastSquares(); range=nothing, sill=nothing, nugget=nothing)

Fit theoretical variogram type `G` to empirical variogram `g` using algorithm `algo`.

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
fit(G::Type{<:Variogram}, g::EmpiricalVariogram, algo::FitAlgo=WeightedLeastSquares(); kwargs...) =
  fit_impl(G, g, algo; kwargs...) |> first

"""
    fit(Gs, g, algo=WeightedLeastSquares(); kwargs...)

Fit theoretical variogram types `Gs` to empirical variogram `g`
using algorithm `algo` and return the one with minimum error.

## Examples

```julia
julia> fit([SphericalVariogram, ExponentialVariogram], g)
```
"""
function fit(Gs, g::EmpiricalVariogram, algo::FitAlgo=WeightedLeastSquares(); kwargs...)
  # fit each variogram type
  res = [fit_impl(G, g, algo; kwargs...) for G in Gs]
  γs, ϵs = first.(res), last.(res)

  # return best candidate
  γs[argmin(ϵs)]
end

"""
    fit(Variogram, g, algo=WeightedLeastSquares(); kwargs...)

Fit all stationary `Variogram` models to empirical variogram `g`
using algorithm `algo` and return the one with minimum error.

## Examples

```julia
julia> fit(Variogram, g)
julia> fit(Variogram, g, WeightedLeastSquares())
```
"""
fit(::Type{Variogram}, g::EmpiricalVariogram, algo::FitAlgo=WeightedLeastSquares(); kwargs...) =
  fit(stationaryvariograms(), g, algo; kwargs...)

"""
    fit(G, g, weightfun; kwargs...)
    fit(Variogram, g, weightfun; kwargs...)

Convenience method that forwards the weighting function `weightfun`
to the `WeightedLeastSquares` algorithm.

## Examples

```julia
fit(SphericalVariogram, g, h -> exp(-h))
fit(Variogram, g, h -> exp(-h/100))
```
"""
fit(G, g::EmpiricalVariogram, weightfun::Function; kwargs...) = fit(G, g, WeightedLeastSquares(weightfun); kwargs...)

# ---------------
# IMPLEMENTATION
# ---------------

function fit_impl(
  G::Type{<:Variogram},
  g::EmpiricalVariogram,
  algo::WeightedLeastSquares;
  range=nothing,
  sill=nothing,
  nugget=nothing,
  maxrange=nothing,
  maxsill=nothing,
  maxnugget=nothing
)
  # custom ball of given radius
  ball(r) = MetricBall(r, g.distance)

  # coordinates of empirical variogram
  x = g.abscissas
  y = g.ordinates
  n = g.counts

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
  w = _weights(algo.weightfun, x, n)

  # objective function
  function J(θ)
    γ = G(ball(θ[1]), sill=θ[2], nugget=θ[3])
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
  nₒ = isnothing(nugget′) ? 0.01 * smax : nugget′
  θₒ = [rₒ, sₒ, nₒ]

  # box constraints
  δ = 1e-8
  rₗ, rᵤ = isnothing(range′) ? (zero(rmax), rmax) : (range′ - δ, range′ + δ)
  sₗ, sᵤ = isnothing(sill′) ? (zero(smax), smax) : (sill′ - δ, sill′ + δ)
  nₗ, nᵤ = isnothing(nugget′) ? (zero(nmax), nmax) : (nugget′ - δ, nugget′ + δ)
  l = [rₗ, sₗ, nₗ]
  u = [rᵤ, sᵤ, nᵤ]

  # solve optimization problem
  θ, ϵ = _optimize(J, L, λ, l, u, θₒ)

  # optimal variogram (with units)
  γ = G(ball(θ[1] * ux), sill=θ[2] * uy, nugget=θ[3] * uy)

  γ, ϵ
end

function fit_impl(
  G::Type{<:PowerVariogram},
  g::EmpiricalVariogram,
  algo::WeightedLeastSquares;
  scaling=nothing,
  nugget=nothing,
  exponent=nothing,
  maxscaling=nothing,
  maxnugget=nothing,
  maxexponent=nothing
)
  # coordinates of empirical variogram
  x = g.abscissas
  y = g.ordinates
  n = g.counts

  # discard invalid bins
  x = x[n .> 0]
  y = y[n .> 0]
  n = n[n .> 0]

  # strip units of coordinates
  uy = unit(eltype(y))
  x′ = ustrip.(x)
  y′ = ustrip.(y)

  # strip units of kwargs
  scaling′ = isnothing(scaling) ? scaling : ustrip(uy, scaling)
  nugget′ = isnothing(nugget) ? nugget : ustrip(uy, nugget)
  exponent′ = exponent
  maxscaling′ = isnothing(maxscaling) ? maxscaling : ustrip(uy, maxscaling)
  maxnugget′ = isnothing(maxnugget) ? maxnugget : ustrip(uy, maxnugget)
  maxexponent′ = maxexponent

  # evaluate weights
  w = _weights(algo.weightfun, x, n)

  # objective function
  function J(θ)
    γ = G(scaling=θ[1], nugget=θ[2], exponent=θ[3])
    sum(i -> w[i] * (γ(x′[i]) - y′[i])^2, eachindex(w, x′, y′))
  end

  # linear constraints
  # 1. scaling ≥ 0
  # 2. 0 ≤ exponent ≤ 2
  L(θ) = θ[1] ≥ 0.0 ? 0.0 : -θ[1] + θ[3] ≥ 0.0 ? 0.0 : -θ[3] + 2.0 ≥ θ[3] ? 0.0 : θ[3] - 2.0

  # penalty for linear constraint (J + λL)
  λ = sum(yᵢ -> yᵢ^2, y′)

  # maximum scaling, nugget and exponent
  ymax = maximum(y′)
  smax = isnothing(maxscaling′) ? ymax : maxscaling′
  nmax = isnothing(maxnugget′) ? ymax : maxnugget′
  emax = isnothing(maxexponent′) ? 2.0 : maxexponent′

  # initial guess
  sₒ = isnothing(scaling′) ? smax / 3 : scaling′
  nₒ = isnothing(nugget′) ? 0.01 * nmax : nugget′
  eₒ = isnothing(exponent′) ? 0.95 * emax : exponent′
  θₒ = [sₒ, nₒ, eₒ]

  # box constraints
  δ = 1e-8
  sₗ, sᵤ = isnothing(scaling′) ? (zero(smax), smax) : (scaling′ - δ, scaling′ + δ)
  nₗ, nᵤ = isnothing(nugget′) ? (zero(nmax), nmax) : (nugget′ - δ, nugget′ + δ)
  eₗ, eᵤ = isnothing(exponent′) ? (zero(emax), emax) : (exponent′ - δ, exponent′ + δ)
  l = [sₗ, nₗ, eₗ]
  u = [sᵤ, nᵤ, eᵤ]

  # solve optimization problem
  θ, ϵ = _optimize(J, L, λ, l, u, θₒ)

  # optimal variogram (with units)
  γ = G(scaling=θ[1] * uy, nugget=θ[2] * uy, exponent=θ[3])

  γ, ϵ
end

_weights(f, x, n) = isnothing(f) ? n / sum(n) : map(xᵢ -> ustrip(f(xᵢ)), x)

function _optimize(J, L, λ, l, u, θₒ)
  sol = Optim.optimize(θ -> J(θ) + λ * L(θ), l, u, θₒ)
  ϵ = Optim.minimum(sol)
  θ = Optim.minimizer(sol)
  θ, ϵ
end
