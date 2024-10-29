# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PiecewiseLinearTransiogram(abscissas, ordinates)
    PiecewiseLinearTransiogram(ball, abscissas, ordinates)

A piecewise-linear transiogram model with `abscissas` and
matrix `ordinates` obtained from an [`EmpiricalTransiogram`](@ref).
Optionally, specify a metric `ball` to model anisotropy.

## References

* Li, W. & Zhang, C. 2005. [Application of Transiograms to Markov Chain
  Simulation and Spatial Uncertainty Assessment of Land-Cover Classes]
  (https://www.tandfonline.com/doi/abs/10.2747/1548-1603.42.4.297)

* Li, W. & Zhang, C. 2010. [Linear interpolation and joint model fitting
  of experimental transiograms for Markov chain simulation of categorical
  spatial variables](https://www.tandfonline.com/doi/abs/10.1080/13658810903127991)
"""
struct PiecewiseLinearTransiogram{ℒ<:Len,M,B<:MetricBall} <: Transiogram
  abscissas::Vector{ℒ}
  ordinates::Vector{M}
  ordinfinity::M
  ball::B
end

function PiecewiseLinearTransiogram(ball::MetricBall, abscissas::AbstractVector, ordinates::AbstractMatrix)
  n = length(abscissas)
  m = size(ordinates, 1)

  # abscissa vector
  x = abscissas

  # ordinate matrices
  Y = map(1:n) do k
    SMatrix{m,m}(ordinates[i, j][k] for i in 1:m, j in 1:m)
  end

  # ordinate matrix at infinity
  p = normalize(diag(last(Y)), 1)
  Y∞ = SMatrix{m,m}(p[j] for i in 1:m, j in 1:m)

  # metric ball
  b = ball

  PiecewiseLinearTransiogram(x, Y, Y∞, b)
end

function PiecewiseLinearTransiogram(abscissas::AbstractVector, ordinates::AbstractMatrix)
  ball = MetricBall(oneunit(eltype(abscissas)))
  PiecewiseLinearTransiogram(ball, abscissas, ordinates)
end

function (t::PiecewiseLinearTransiogram)(h)
  hs = t.abscissas
  if h < first(hs) # left extrapolation
    ((first(hs) - h) * I + h * first(t.ordinates)) / first(hs)
  elseif h > last(hs) # right extrapolation
    t.ordinfinity
  else # middle interpolation
    k = findlast(<(h), hs)
    ((hs[k + 1] - h) * t.ordinates[k] + (h - hs[k]) * t.ordinates[k + 1]) / (hs[k + 1] - hs[k])
  end
end
