# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MatrixExponentialTransiogram(rate)

An exponential transiogram with transition `rate` matrix.

    MatrixExponentialTransiogram(lengths, proportions)

Alternatively, build transition rate matrix from mean `lengths`
and relative `proportions`.

    MatrixExponentialTransiogram(ball, rate)

Alternatively, use a custom metric `ball`.

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)

* Carle et al 1998. [Conditional Simulation of Hydrofacies Architecture:
  A Transition Probability/Markov Approach](https://doi.org/10.2110/sepmcheg.01.147)
"""
struct MatrixExponentialTransiogram{B<:MetricBall,R<:StaticMatrix} <: Transiogram
  ball::B
  rate::R

  function MatrixExponentialTransiogram{B,R}(ball, rate) where {B<:MetricBall,R<:StaticMatrix}
    if !allequal(size(rate))
      throw(ArgumentError("transition rate matrix must be square"))
    end
    new(ball, rate)
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

MatrixExponentialTransiogram(ball::MetricBall, lens::Tuple, props::Tuple) =
  MatrixExponentialTransiogram(ball, baseratematrix(lens, props))

MatrixExponentialTransiogram(lens::Tuple, props::Tuple) = MatrixExponentialTransiogram(baseratematrix(lens, props))

MatrixExponentialTransiogram() = MatrixExponentialTransiogram((1.0u"m", 1.0u"m"), (0.5, 0.5))

meanlengths(t::MatrixExponentialTransiogram) = Tuple(1 ./ -diag(t.rate))

(t::MatrixExponentialTransiogram)(h) = exp(aslen(h) * t.rate)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, t::MatrixExponentialTransiogram)
  ioctx = IOContext(io, :compact => true)
  summary(ioctx, t)
  print(ioctx, "(")
  print(ioctx, "rate: ")
  _printvec(ioctx, t.rate, 1)
  if isisotropic(t.ball)
    print(ioctx, ", range: ", first(radii(t.ball)))
  else
    print(ioctx, ", ranges: ", radii(t.ball))
  end
  m = nameof(typeof(metric(t.ball)))
  print(ioctx, ", distance: ", m)
  print(ioctx, ")")
end

function Base.show(io::IO, ::MIME"text/plain", t::MatrixExponentialTransiogram)
  ioctx = IOContext(io, :compact => true, :limit => true)
  summary(ioctx, t)
  println(ioctx)
  print(ioctx, "├─ rate: ")
  _printlnvec(ioctx, t.rate, 3)
  if isisotropic(t.ball)
    println(ioctx, "├─ range: ", first(radii(t.ball)))
  else
    println(ioctx, "├─ ranges: ", radii(t.ball))
  end
  m = nameof(typeof(metric(t.ball)))
  print(ioctx, "└─ distance: ", m)
end

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
  SMatrix{nₗ,nₗ}(i == j ? -1 / l[i] : (p[j] / (1 - p[i])) / l[i] for i in 1:nₗ, j in 1:nₗ)
end
