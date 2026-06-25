# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _fit(
  T::Type{<:Transiogram},
  t::EmpiricalTransiogram,
  algo::WeightedLeastSquares;
  range=nothing,
  proportions=nothing,
  maxrange=nothing
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

  # number of categories and free logits
  k = size(Y, 1)
  m = k - 1

  # objective function (θ = [range, logits...])
  function J(θ)
    b = ball(θ[1])
    p = isnothing(proportions) ? _softmax(θ[2:end]) : proportions
    τ = T(b, proportions=p)
    mat(i) = getindex.(Y, i)
    err(i) = sum(abs2, τ(x′[i]) - mat(i))
    sum(i -> w[i] * err(i), eachindex(w, x′))
  end

  # maximum range
  rmax = isnothing(maxrange′) ? maximum(x′) : maxrange′

  # initial guess and box constraints for the range
  δ = oftype(rmax, 1e-8)
  rₒ = isnothing(range′) ? rmax / 3 : range′
  rₗ, rᵤ = isnothing(range′) ? (δ, rmax) : (range′ - δ, range′ + δ)

  # logits are unconstrained, proportions are recovered via softmax
  λₒ = ntuple(i -> oftype(rₒ, 0.0), m)
  λₗ = ntuple(i -> oftype(rₗ, -Inf), m)
  λᵤ = ntuple(i -> oftype(rᵤ, Inf), m)

  θₒ = isnothing(proportions) ? [rₒ, λₒ...] : [rₒ]
  l = isnothing(proportions) ? [rₗ, λₗ...] : [rₗ]
  u = isnothing(proportions) ? [rᵤ, λᵤ...] : [rᵤ]

  # solve optimization problem
  θ, ϵ = _optimize(J, l, u, θₒ)

  # optimal transiogram (with units)
  bopt = ball(θ[1] * ux)
  popt = isnothing(proportions) ? _softmax(θ[2:end]) : proportions
  τ = T(bopt, proportions=popt)

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
  maxlengths=nothing
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

  # number of categories and free logits
  k = size(Y, 1)
  m = k - 1

  # objective function (θ = [range, lengths..., logits...])
  function J(θ)
    b = ball(θ[1])
    l = _lengths(θ[2:(k + 1)])
    p = isnothing(proportions) ? _softmax(θ[(k + 2):end]) : proportions
    τ = T(b, lengths=l, proportions=p)
    mat(i) = getindex.(Y, i)
    err(i) = sum(abs2, τ(x′[i]) - mat(i))
    sum(i -> w[i] * err(i), eachindex(w, x′))
  end

  # maximum range and lengths
  xmax = maximum(x′)
  rmax = isnothing(maxrange′) ? xmax : maxrange′
  lmax = isnothing(maxlengths′) ? ntuple(i -> xmax, k) : maxlengths′

  # initial guess and box constraints for the range and lengths
  δ = oftype(rmax, 1e-8)
  rₒ = isnothing(range′) ? rmax / 3 : range′
  lₒ = isnothing(lengths′) ? lmax ./ 3 : lengths′
  rₗ, rᵤ = isnothing(range′) ? (δ, rmax) : (range′ - δ, range′ + δ)
  lₗ, lᵤ = isnothing(lengths′) ? (ntuple(i -> δ, k), lmax) : (lengths′ .- δ, lengths′ .+ δ)

  # logits are unconstrained, proportions are recovered via softmax
  λₗ = ntuple(i -> oftype(rₗ, -Inf), m)
  λᵤ = ntuple(i -> oftype(rᵤ, Inf), m)

  # collect all box constraints for range, lengths and logits
  l = isnothing(proportions) ? [rₗ, lₗ..., λₗ...] : [rₗ, lₗ...]
  u = isnothing(proportions) ? [rᵤ, lᵤ..., λᵤ...] : [rᵤ, lᵤ...]

  # the matrix-exponential misfit is non-convex, so a single start can land in
  # a poor local minimum. Run a small deterministic multistart that varies the
  # free blocks (length scale and proportion logits) and keep the best fit.
  lstarts = isnothing(lengths′) ? [s .* lmax for s in (0.1, 0.25, 0.5, 0.75, 0.9)] : [lₒ]
  λstarts = isnothing(proportions) ? [ntuple(c -> b * (c - m / 2), m) for b in (-1.0, 0.0, 1.0)] : [()]
  starts = [[rₒ, lₒ..., λₒ...] for lₒ in lstarts for λₒ in λstarts]
  sols = [_optimize(J, l, u, θₒ) for θₒ in starts]
  θ, ϵ = sols[argmin(last.(sols))]

  # optimal transiogram (with units)
  bopt = ball(θ[1] * ux)
  lopt = _lengths(θ[2:(k + 1)] * ux)
  popt = isnothing(proportions) ? _softmax(θ[(k + 2):end]) : proportions
  τ = T(bopt, lengths=lopt, proportions=popt)

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

function _softmax(θ)
  T = eltype(θ)
  n = length(θ)
  e = ntuple(i -> i ≤ n ? exp(θ[i]) : one(T), n + 1)
  s = sum(e)
  ntuple(i -> e[i] / s, n + 1)
end

_lengths(θ) = ntuple(i -> θ[i], length(θ))
