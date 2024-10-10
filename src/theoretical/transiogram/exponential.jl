# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExponentialTransiogram(rate; levels=l)
    ExponentialTransiogram(ball, rate; levels=l)

An exponential transiogram with transition `rate` matrix.
Optionally, specify a metric `ball` to model anisotropy,
and the `levels` or categories.

    ExponentialTransiogram(lengths, proportions; levels=l)
    ExponentialTransiogram(ball, lengths, proportions; levels=l)

Alternatively, build transition rate matrix from mean `lengths`
and relative `proportions`.

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)

* Carle et al 1998. [Conditional Simulation of Hydrofacies Architecture:
  A Transition Probability/Markov Approach](https://doi.org/10.2110/sepmcheg.01.147)
"""
struct ExponentialTransiogram{R<:StaticMatrix,L<:AbstractVector,B<:MetricBall} <: Transiogram
  rate::R
  levs::L
  ball::B

  function ExponentialTransiogram{R,L,B}(rate, levs, ball) where {R<:StaticMatrix,L<:AbstractVector,B<:MetricBall}
    if !allequal(size(rate))
      throw(ArgumentError("transition rate matrix must be square"))
    end
    if length(levs) != size(rate, 1)
      throw(ArgumentError("levels do not match size of transition rate matrix"))
    end
    new(rate, levs, ball)
  end
end

function ExponentialTransiogram(ball::MetricBall, rate::AbstractMatrix; levels=1:size(rate, 1))
  srate = SMatrix{size(rate)...}(rate)
  ExponentialTransiogram{typeof(srate),typeof(levels),typeof(ball)}(srate, levels, ball)
end

function ExponentialTransiogram(rate::AbstractMatrix; levels=1:size(rate, 1))
  ball = MetricBall(1 / unit(eltype(rate)))
  ExponentialTransiogram(ball, rate; levels)
end

ExponentialTransiogram(ball::MetricBall, lens::AbstractVector, props::AbstractVector; levels=1:length(lens)) =
  ExponentialTransiogram(ball, baseratematrix(lens, props); levels)

ExponentialTransiogram(lens::AbstractVector, props::AbstractVector; levels=1:length(lens)) =
  ExponentialTransiogram(baseratematrix(lens, props); levels)

ranges(t::Transiogram) = 1 ./ -diag(t.rate)

levels(t::ExponentialTransiogram) = t.levs

(t::ExponentialTransiogram)(h) = exp(h * t.rate)

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
