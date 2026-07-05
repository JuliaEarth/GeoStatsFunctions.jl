# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

_fit(G::Type{<:Variogram}, g::EmpiricalVariogram, w::Union{Function,Nothing}; kwargs...) =
  nvariables(g) == 1 ? _fitunivariate(G, g, w; kwargs...) : _fitmultivariate(G, g, w; kwargs...)

function _fitunivariate(
  G::Type{<:Variogram},
  g::EmpiricalVariogram,
  w::Union{Function,Nothing};
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
  y = g.ordinates[1]
  n = g.counts

  # discard invalid bins
  x = x[n .> 0]
  y = y[n .> 0]
  n = n[n .> 0]

  # evaluate weights
  ω = _weights(w, x, n)

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

  # maximum range, sill and nugget
  xmax = maximum(x′)
  ymax = maximum(y′)
  rmax = isnothing(maxrange′) ? xmax : maxrange′
  smax = isnothing(maxsill′) ? ymax : maxsill′
  nmax = isnothing(maxnugget′) ? ymax : maxnugget′

  # objective function
  function J(θ)
    γ = G(ball(θ[1]), sill=θ[2], nugget=θ[3])
    sum(i -> ω[i] * (γ(x′[i]) - y′[i])^2, eachindex(ω, x′, y′))
  end

  # linear constraint (sill ≥ nugget)
  L(θ) = θ[2] ≥ θ[3] ? 0.0 : θ[3] - θ[2]

  # penalty for linear constraint (J + λL)
  λ = sum(abs2, y′)

  # box constraints
  δ = oftype(rmax, 1e-8)
  rₗ, rᵤ = isnothing(range′) ? (δ, rmax) : (range′, range′)
  sₗ, sᵤ = isnothing(sill′) ? (zero(smax), smax) : (sill′, sill′)
  nₗ, nᵤ = isnothing(nugget′) ? (zero(nmax), nmax) : (nugget′, nugget′)
  θₗ = [rₗ, sₗ, nₗ]
  θᵤ = [rᵤ, sᵤ, nᵤ]

  # initial guess
  rₒ = isnothing(range′) ? rmax / 3 : range′
  sₒ = isnothing(sill′) ? 0.95 * smax : sill′
  nₒ = isnothing(nugget′) ? 0.01 * smax : nugget′
  θₒ = [rₒ, sₒ, nₒ]

  # solve optimization problem
  θ, ϵ = _optimize(θ -> J(θ) + λ * L(θ), θₗ, θᵤ, θₒ)

  # optimal variogram (with units)
  γ = G(ball(θ[1] * ux), sill=θ[2] * uy, nugget=θ[3] * uy)

  γ, ϵ
end

function _fitunivariate(
  G::Type{<:PowerVariogram},
  g::EmpiricalVariogram,
  w::Union{Function,Nothing};
  scaling=nothing,
  nugget=nothing,
  exponent=nothing,
  maxscaling=nothing,
  maxnugget=nothing,
  maxexponent=nothing
)
  # coordinates of empirical variogram
  x = g.abscissas
  y = g.ordinates[1]
  n = g.counts

  # discard invalid bins
  x = x[n .> 0]
  y = y[n .> 0]
  n = n[n .> 0]

  # evaluate weights
  ω = _weights(w, x, n)

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

  # maximum scaling, nugget and exponent
  ymax = maximum(y′)
  smax = isnothing(maxscaling′) ? ymax : maxscaling′
  nmax = isnothing(maxnugget′) ? ymax : maxnugget′
  emax = isnothing(maxexponent′) ? 2.0 : maxexponent′

  # objective function
  function J(θ)
    γ = G(scaling=θ[1], nugget=θ[2], exponent=θ[3])
    sum(i -> ω[i] * (γ(x′[i]) - y′[i])^2, eachindex(ω, x′, y′))
  end

  # linear constraints
  # 1. scaling ≥ 0
  # 2. 0 ≤ exponent ≤ 2
  L(θ) = θ[1] ≥ 0.0 ? 0.0 : -θ[1] + θ[3] ≥ 0.0 ? 0.0 : -θ[3] + 2.0 ≥ θ[3] ? 0.0 : θ[3] - 2.0

  # penalty for linear constraint (J + λL)
  λ = sum(abs2, y′)

  # box constraints
  δ = oftype(smax, 1e-8)
  sₗ, sᵤ = isnothing(scaling′) ? (δ, smax) : (scaling′, scaling′)
  nₗ, nᵤ = isnothing(nugget′) ? (zero(nmax), nmax) : (nugget′, nugget′)
  eₗ, eᵤ = isnothing(exponent′) ? (zero(emax), emax) : (exponent′, exponent′)
  θₗ = [sₗ, nₗ, eₗ]
  θᵤ = [sᵤ, nᵤ, eᵤ]

  # initial guess
  sₒ = isnothing(scaling′) ? smax / 3 : scaling′
  nₒ = isnothing(nugget′) ? 0.01 * nmax : nugget′
  eₒ = isnothing(exponent′) ? 0.95 * emax : exponent′
  θₒ = [sₒ, nₒ, eₒ]

  # solve optimization problem
  θ, ϵ = _optimize(θ -> J(θ) + λ * L(θ), θₗ, θᵤ, θₒ)

  # optimal variogram (with units)
  γ = G(scaling=θ[1] * uy, nugget=θ[2] * uy, exponent=θ[3])

  γ, ϵ
end

function _fitmultivariate(
  G::Type{<:Variogram},
  g::EmpiricalVariogram,
  w::Union{Function,Nothing};
  range=nothing,
  sill=nothing,
  nugget=nothing,
  maxrange=nothing,
  maxsill=nothing,
  maxnugget=nothing
)
  # number of variables
  k = nvariables(g)

  # custom ball of given radius
  ball(r) = MetricBall(r, g.distance)

  # coordinates of empirical variogram
  x = g.abscissas
  Y = g.ordinates
  n = g.counts

  # discard invalid bins
  x = x[n .> 0]
  Y = [y[n .> 0] for y in Y]
  n = n[n .> 0]

  # evaluate weights
  ω = _weights(w, x, n)

  # strip units of coordinates
  ux = unit(eltype(x))
  uY = [unit(eltype(y)) for y in Y]
  x′ = ustrip.(x)
  Y′ = [ustrip.(y) for y in Y]

  # strip units of kwargs
  range′ = isnothing(range) ? range : _ustrip(ux, range)
  sill′ = isnothing(sill) ? sill : _ustrip(uY, sill)
  nugget′ = isnothing(nugget) ? nugget : _ustrip(uY, nugget)
  maxrange′ = isnothing(maxrange) ? maxrange : _ustrip(ux, maxrange)
  maxsill′ = isnothing(maxsill) ? maxsill : _ustrip(uY, maxsill)
  maxnugget′ = isnothing(maxnugget) ? maxnugget : _ustrip(uY, maxnugget)

  # maximum range, sill and nugget
  xmax = maximum(x′)
  Ymax = [maximum(y′) for y′ in Y′]
  rmax = isnothing(maxrange′) ? xmax : maxrange′
  smax = isnothing(maxsill′) ? _posmat(Ymax) : maxsill′
  nmax = isnothing(maxnugget′) ? _posmat(Ymax) : maxnugget′

  # number of parameters in lower triangular matrix
  p = k * (k + 1) ÷ 2

  # objective function
  function J(θ)
    r = θ[1]
    Bₒ = _coefmat(θ[2:(p + 1)], k)
    B₁ = _coefmat(θ[(p + 2):end], k)
    γ = Bₒ * NuggetEffect() + B₁ * G(ball(r))
    sum(eachindex(ω, x′)) do i
      ωᵢ = ω[i]
      Γᵢ = γ(x′[i])
      sum(j -> ωᵢ * (Γᵢ[j] - Y′[j][i])^2, eachindex(Γᵢ, Y′))
    end
  end

  # linear constraint (sill ≥ nugget)
  function L(θ)
    Bₒ = _coefmat(θ[2:(p + 1)], k)
    B₁ = _coefmat(θ[(p + 2):end], k)
    sum(zip(Bₒ, B₁)) do (bₒ, b₁)
      b₁ ≥ bₒ ? 0.0 : bₒ - b₁
    end
  end

  # penalty for linear constraint (J + λL)
  λ = sum(sum(abs2, y′) for y′ in Y′)

  # box constraints
  δ = oftype(rmax, 1e-8)
  rₗ, rᵤ = isnothing(range′) ? (δ, rmax) : (range′, range′)
  nₗ, nᵤ = isnothing(nugget′) ? (_coefvec(δ*I(k)), _coefvec(nmax)) : (_coefvec(nugget′), _coefvec(nugget′))
  sₗ, sᵤ = isnothing(sill′) ? (_coefvec(δ*I(k)), _coefvec(smax)) : (_coefvec(sill′), _coefvec(sill′))
  θₗ = [rₗ, nₗ..., sₗ...]
  θᵤ = [rᵤ, nᵤ..., sᵤ...]

  # initial guess
  rₒ = isnothing(range′) ? rmax / 3 : range′
  nₒ = isnothing(nugget′) ? _coefvec(0.01 * smax) : _coefvec(nugget′)
  sₒ = isnothing(sill′) ? _coefvec(0.95 * smax) : _coefvec(sill′)
  θₒ = [rₒ, nₒ..., sₒ...]

  # solve optimization problem
  θ, ϵ = _optimize(θ -> J(θ) + λ * L(θ), θₗ, θᵤ, θₒ)

  # optimal variogram (with units)
  γ = let
    r = θ[1] * ux
    Bₒ = _coefmat(θ[2:(p + 1)], k) .* uY
    B₁ = _coefmat(θ[(p + 2):end], k) .* uY
    Bₒ * NuggetEffect() + B₁ * G(ball(r))
  end

  γ, ϵ
end
