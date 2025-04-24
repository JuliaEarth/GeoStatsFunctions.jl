# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MatrixExponentialTransiogram(ball, rate)

A matrix-exponential transiogram with given metric `ball` and
transition `rate` matrix.

    MatrixExponentialTransiogram(ball; lengths, proportions)

Alternatively, build transition `rate` matrix from nominal `lengths`
and relative `proportions`.

    MatrixExponentialTransiogram(rate; range)
    MatrixExponentialTransiogram(rate; ranges, rotation)

Alternatively, build metric `ball` from single `range` or from
multiple `ranges` and `rotation` matrix.

    MatrixExponentialTransiogram(; range, lengths, proportions)
    MatrixExponentialTransiogram(; ranges, rotation, lengths, proportions)

Alternatively, build metric `ball` and transition `rate` matrix
from combination of previous parameters.

## Examples

```julia
MatrixExponentialTransiogram(lengths=(3.0, 2.0, 1.0), proportions=(1/3, 1/3, 1/3))
MatrixExponentialTransiogram(ranges=(2.0, 1.0))
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

MatrixExponentialTransiogram(ball::MetricBall; lengths=(1.0, 1.0), proportions=(0.5, 0.5)) =
  MatrixExponentialTransiogram(ball, baseratematrix(lengths, proportions))

MatrixExponentialTransiogram(rate::AbstractMatrix; range=1.0, ranges=nothing, rotation=I) =
  MatrixExponentialTransiogram(rangeball(range, ranges, rotation), rate)

MatrixExponentialTransiogram(; range=1.0, ranges=nothing, rotation=I, lengths=(1.0, 1.0), proportions=(0.5, 0.5)) =
  MatrixExponentialTransiogram(rangeball(range, ranges, rotation), baseratematrix(lengths, proportions))

constructor(::MatrixExponentialTransiogram) = MatrixExponentialTransiogram

Base.range(t::MatrixExponentialTransiogram) = maximum(meanlengths(t))

scale(t::MatrixExponentialTransiogram, s::Real) = MatrixExponentialTransiogram(s * t.ball, t.rate)

meanlengths(t::MatrixExponentialTransiogram) = Tuple(1 ./ -diag(t.rate)) .* ustrip(maximum(radii(t.ball)))

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
