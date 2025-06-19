# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _fit(
  T::Type{<:Transiogram},
  t::EmpiricalTransiogram,
  algo::WeightedLeastSquares;
  range=nothing,
  proportions=nothing,
  maxrange=nothing,
  maxproportions=nothing
)
  # custom ball of given radius
  ball(r) = MetricBall(r, t.distance)

  # coordinates of empirical transiogram
  x = t.abscissas
  Y = t.ordinates
  n = t.counts

  # discard invalid bins
  x = x[n .> 0]
  Y = [y[n .> 0] for y in Y]
  n = n[n .> 0]

  # strip units of coordinates
  ux = unit(eltype(x))
  x′ = ustrip.(x)

  # strip units of kwargs
  range′ = isnothing(range) ? range : _ustrip(ux, range)
  maxrange′ = isnothing(maxrange) ? maxrange : _ustrip(ux, maxrange)

  # evaluate weights
  w = _weights(algo.weightfun, x, n)

  # auxiliary variables
  k = size(Y, 1)
  V = eltype(first(Y))
  _ones = ntuple(i -> one(V), k)
  _zeros = ntuple(i -> zero(V), k)

  # objective function
  function J(θ)
    b = ball(θ[1])
    p = _proportions(θ[2:end])
    τ = T(b, proportions=p)
    mat(i) = getindex.(Y, i)
    err(i) = sum(abs2, τ(x′[i]) - mat(i))
    sum(i -> w[i] * err(i), eachindex(w, x′))
  end

  # linear constraint (sum(proportions) == 1)
  L(θ) = abs(sum(θ[2:end]) - 1)

  # penalty for linear constraint (J + λL)
  λ = sum(y -> sum(abs2, y), Y)

  # maximum range and proportions
  xmax = maximum(x′)
  ymax = _ones
  rmax = isnothing(maxrange′) ? xmax : maxrange′
  pmax = isnothing(maxproportions) ? ymax : maxproportions

  # initial guess
  rₒ = isnothing(range′) ? rmax / 3 : range′
  pₒ = isnothing(proportions) ? 0.95 .* pmax : proportions
  θₒ = [rₒ, pₒ...]

  # box constraints
  δ = 1e-8
  rₗ, rᵤ = isnothing(range′) ? (zero(rmax), rmax) : (range′ - δ, range′ + δ)
  pₗ, pᵤ = isnothing(proportions) ? (_zeros, pmax) : (proportions .- δ, proportions .+ δ)
  l = [rₗ, pₗ...]
  u = [rᵤ, pᵤ...]

  # solve optimization problem
  θ, ϵ = _optimize(J, L, λ, l, u, θₒ)

  # optimal transiogram (with units)
  b = ball(θ[1] * ux)
  p = _proportions(θ[2:end])
  τ = T(b, proportions=p)

  τ, ϵ
end

function _fit(
  T::Type{<:MatrixExponentialTransiogram},
  t::EmpiricalTransiogram,
  algo::WeightedLeastSquares;
  range=nothing,
  lengths=nothing,
  proportions=nothing,
  maxrange=nothing,
  maxlengths=nothing,
  maxproportions=nothing
)
  # custom ball of given radius
  ball(r) = MetricBall(r, t.distance)

  # coordinates of empirical transiogram
  x = t.abscissas
  Y = t.ordinates
  n = t.counts

  # discard invalid bins
  x = x[n .> 0]
  Y = [y[n .> 0] for y in Y]
  n = n[n .> 0]

  # strip units of coordinates
  ux = unit(eltype(x))
  x′ = ustrip.(x)

  # strip units of kwargs
  range′ = isnothing(range) ? range : _ustrip(ux, range)
  lengths′ = isnothing(lengths) ? lengths : _ustrip.(ux, lengths)
  maxrange′ = isnothing(maxrange) ? maxrange : _ustrip(ux, maxrange)
  maxlengths′ = isnothing(maxlengths) ? maxlengths : _ustrip.(ux, maxlengths)

  # evaluate weights
  w = _weights(algo.weightfun, x, n)

  # auxiliary variables
  k = size(Y, 1)
  V = eltype(first(Y))
  _ones = ntuple(i -> one(V), k)
  _zeros = ntuple(i -> zero(V), k)

  # objective function
  function J(θ)
    b = ball(θ[1])
    l = _lengths(θ[2:(k + 1)])
    p = _proportions(θ[(k + 2):end])
    τ = T(b, lengths=l, proportions=p)
    mat(i) = getindex.(Y, i)
    err(i) = sum(abs2, τ(x′[i]) - mat(i))
    sum(i -> w[i] * err(i), eachindex(w, x′))
  end

  # linear constraint (sum(proportions) == 1)
  L(θ) = abs(sum(θ[(k + 2):end]) - 1)

  # penalty for linear constraint (J + λL)
  λ = sum(y -> sum(abs2, y), Y)

  # maximum range, proportions and lengths
  xmax = maximum(x′)
  ymax = _ones
  rmax = isnothing(maxrange′) ? xmax : maxrange′
  lmax = isnothing(maxlengths′) ? ntuple(i -> xmax, k) : maxlengths′
  pmax = isnothing(maxproportions) ? ymax : maxproportions

  # initial guess
  rₒ = isnothing(range′) ? rmax / 3 : range′
  lₒ = isnothing(lengths′) ? lmax ./ 3 : lengths′
  pₒ = isnothing(proportions) ? 0.95 .* pmax : proportions
  θₒ = [rₒ, lₒ..., pₒ...]

  # box constraints
  δ = 1e-8
  rₗ, rᵤ = isnothing(range′) ? (zero(rmax), rmax) : (range′ - δ, range′ + δ)
  lₗ, lᵤ = isnothing(lengths′) ? (_zeros, lmax) : (lengths′ .- δ, lengths′ .+ δ)
  pₗ, pᵤ = isnothing(proportions) ? (_zeros, pmax) : (proportions .- δ, proportions .+ δ)
  l = [rₗ, lₗ..., pₗ...]
  u = [rᵤ, lᵤ..., pᵤ...]

  # solve optimization problem
  θ, ϵ = _optimize(J, L, λ, l, u, θₒ)

  # optimal transiogram (with units)
  b = ball(θ[1] * ux)
  l = _lengths(θ[2:(k + 1)] * ux)
  p = _proportions(θ[(k + 2):end])
  τ = T(b, lengths=l, proportions=p)

  τ, ϵ
end

function _fit(T::Type{<:PiecewiseLinearTransiogram}, t::EmpiricalTransiogram, ::WeightedLeastSquares)
  # auxiliary variables
  n = length(t.abscissas)
  k = size(t.ordinates, 1)

  # abscissa vector
  x = t.abscissas

  # ordinate matrices
  Y = map(1:n) do m
    SMatrix{k,k}(t.ordinates[i, j][m] for i in 1:k, j in 1:k)
  end

  # theoretical model
  τ = T(x, Y)

  # zero error
  ϵ = 0.0

  τ, ϵ
end

# -----------------
# HELPER FUNCTIONS
# -----------------

function _proportions(θ)
  p = clamp.(θ, 0, 1)
  s = sum(abs, p)
  ntuple(i -> p[i] / s, length(p))
end

_lengths(θ) = ntuple(i -> θ[i], length(θ))
