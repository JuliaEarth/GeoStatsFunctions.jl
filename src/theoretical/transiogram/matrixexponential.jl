# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MatrixExponentialTransiogram(rate)
    MatrixExponentialTransiogram(ball, rate)

An exponential transiogram with transition `rate` matrix.
Optionally, specify a metric `ball` to model anisotropy.

    MatrixExponentialTransiogram(lengths, proportions)
    MatrixExponentialTransiogram(ball, lengths, proportions)

Alternatively, build transition rate matrix from mean `lengths`
and relative `proportions`.

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)

* Carle et al 1998. [Conditional Simulation of Hydrofacies Architecture:
  A Transition Probability/Markov Approach](https://doi.org/10.2110/sepmcheg.01.147)
"""
struct MatrixExponentialTransiogram{R<:StaticMatrix,B<:MetricBall} <: Transiogram
  rate::R
  ball::B

  function MatrixExponentialTransiogram{R,B}(rate, ball) where {R<:StaticMatrix,B<:MetricBall}
    if !allequal(size(rate))
      throw(ArgumentError("transition rate matrix must be square"))
    end
    new(rate, ball)
  end
end

function MatrixExponentialTransiogram(ball::MetricBall, rate::AbstractMatrix)
  srate = SMatrix{size(rate)...}(float(asinvlen.(rate)))
  MatrixExponentialTransiogram{typeof(srate),typeof(ball)}(srate, ball)
end

function MatrixExponentialTransiogram(rate::AbstractMatrix)
  ball = MetricBall(1 / unit(eltype(rate)))
  MatrixExponentialTransiogram(ball, rate)
end

MatrixExponentialTransiogram(ball::MetricBall, props::AbstractVector) =
  MatrixExponentialTransiogram(ball, baseratematrix(lens, props))

MatrixExponentialTransiogram(lens::AbstractVector, props::AbstractVector) =
  MatrixExponentialTransiogram(baseratematrix(lens, props))

ranges(t::Transiogram) = 1 ./ -diag(t.rate)

(t::MatrixExponentialTransiogram)(h) = exp(h * t.rate)

# -----------------
# HELPER FUNCTIONS
# -----------------

function baseratematrix(l, p)
  nₗ = length(l)
  nₚ = length(p)

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
  map(Iterators.product(1:nₗ, 1:nₗ)) do (i, j)
    if i == j
      -1 / l[i]
    else
      (p[j] / (1 - p[i])) / l[i]
    end
  end
end
