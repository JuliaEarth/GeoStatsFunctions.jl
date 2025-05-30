# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _fit(
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
  range′ = isnothing(range) ? range : _ustrip(ux, range)
  sill′ = isnothing(sill) ? sill : _ustrip(uy, sill)
  nugget′ = isnothing(nugget) ? nugget : _ustrip(uy, nugget)
  maxrange′ = isnothing(maxrange) ? maxrange : _ustrip(ux, maxrange)
  maxsill′ = isnothing(maxsill) ? maxsill : _ustrip(uy, maxsill)
  maxnugget′ = isnothing(maxnugget) ? maxnugget : _ustrip(uy, maxnugget)

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
  λ = sum(abs2, y′)

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

function _fit(
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
  scaling′ = isnothing(scaling) ? scaling : _ustrip(uy, scaling)
  nugget′ = isnothing(nugget) ? nugget : _ustrip(uy, nugget)
  exponent′ = exponent
  maxscaling′ = isnothing(maxscaling) ? maxscaling : _ustrip(uy, maxscaling)
  maxnugget′ = isnothing(maxnugget) ? maxnugget : _ustrip(uy, maxnugget)
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
  λ = sum(abs2, y′)

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
