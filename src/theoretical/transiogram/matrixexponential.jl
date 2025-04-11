# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MatrixExponentialTransiogram(rate)

A matrix-exponential transiogram with transition `rate` matrix.

    MatrixExponentialTransiogram(; range, lengths, proportions)

Alternatively, build Euclidean ball with given `range` and
base transition rate matrix from mean `lengths` and relative
`proportions`.

    MatrixExponentialTransiogram(; ranges, rotation, lengths, proportions)

Alternatively, use multiple `ranges` and `rotation` matrix
to construct an anisotropic model.

    MatrixExponentialTransiogram(ball, rate)

Alternatively, use a custom metric `ball`.

## Examples

```julia
MatrixExponentialTransiogram(lengths=(3.0, 2.0, 1.0))
MatrixExponentialTransiogram(proportions=(0.5, 0.5))
```

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)

* Carle et al 1998. [Conditional Simulation of Hydrofacies Architecture:
  A Transition Probability/Markov Approach](https://doi.org/10.2110/sepmcheg.01.147)

### Notes

The final proportions of the matrix-exponential model can be different
from the relative proportions specified during construction. They are
also a function of the mean lengths.
"""
struct MatrixExponentialTransiogram{B<:MetricBall,R<:StaticMatrix} <: Transiogram
  ball::B
  rate::R

  function MatrixExponentialTransiogram{B,R}(ball, rate) where {B<:MetricBall,R<:StaticMatrix}
    if !allequal(size(rate))
      throw(ArgumentError("transition rate matrix must be square"))
    end
    urate = uconvert.(inv(unit(radius(ball))), rate)
    new{B,typeof(urate)}(ball, urate)
  end
end

function MatrixExponentialTransiogram(ball::MetricBall, rate::AbstractMatrix)
  srate = SMatrix{size(rate)...}(float(asinvlen.(rate)))
  MatrixExponentialTransiogram{typeof(ball),typeof(srate)}(ball, srate)
end

function MatrixExponentialTransiogram(rate::AbstractMatrix)
  ball = MetricBall(1 / unit(eltype(rate)))
  MatrixExponentialTransiogram(ball, rate)
end

MatrixExponentialTransiogram(ball::MetricBall; lengths=(1.0, 1.0), proportions=(0.5, 0.5)) =
  MatrixExponentialTransiogram(ball, baseratematrix(lengths, proportions))

MatrixExponentialTransiogram(; range=1.0, ranges=nothing, rotation=I, lengths=(1.0, 1.0), proportions=(0.5, 0.5)) =
  MatrixExponentialTransiogram(rangeball(range, ranges, rotation), baseratematrix(lengths, proportions))

constructor(::MatrixExponentialTransiogram) = MatrixExponentialTransiogram

function scale(t::MatrixExponentialTransiogram, s::Real)
  T = constructor(t)
  T(s * metricball(t), t.rate)
end

meanlengths(t::MatrixExponentialTransiogram) = Tuple(1 ./ -diag(t.rate))

proportions(t::MatrixExponentialTransiogram) = Tuple(normalize(diag(t(100range(t))), 1))

function (t::MatrixExponentialTransiogram)(h)
  r = radius(t.ball)
  R = ustrip.(t.rate)
  h′, r′ = unitless(h, r)
  v = h′ / r′
  exp(v * R)
end

# -----------------
# HELPER FUNCTIONS
# -----------------

function baseratematrix(l, p)
  nₗ = length(l)
  nₚ = length(p)

  # add units if necessary
  lᵤ = aslen.(l)

  # sanity checks
  if nₗ != nₚ
    throw(ArgumentError("lengths and proportions must have the same length"))
  end
  if !all(pᵢ -> 0 ≤ pᵢ ≤ 1, p)
    throw(ArgumentError("proportions must be in interval [0, 1]"))
  end
  if !(sum(p) ≈ 1)
    throw(ArgumentError("proportions must add up to unit"))
  end

  # Eq. 17 and 18 of Carle et al 1998.
  SMatrix{nₗ,nₗ}(i == j ? -1 / lᵤ[i] : (p[j] / (1 - p[i])) / lᵤ[i] for i in 1:nₗ, j in 1:nₗ)
end
