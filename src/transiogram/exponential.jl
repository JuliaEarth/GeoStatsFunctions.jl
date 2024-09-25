# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExponentialTransiogram(rate, [ball])

An exponential transiogram with transition `rate` matrix.
Optionally, specify a metric `ball` to model anisotropy.

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)
"""
struct ExponentialTransiogram{R<:StaticMatrix,B<:MetricBall} <: Transiogram
  rate::R
  ball::B

  function ExponentialTransiogram{R,B}(rate, ball) where {R<:StaticMatrix,B<:MetricBall}
    if !allequal(size(rate))
      throw(ArgumentError("transition rate matrix must be square"))
    end
    new(rate, ball)
  end
end

function ExponentialTransiogram(rate::AbstractMatrix, ball::MetricBall)
  srate = SMatrix{size(rate)...}(rate)
  ExponentialTransiogram{typeof(srate),typeof(ball)}(srate, ball)
end

function ExponentialTransiogram(rate::AbstractMatrix)
  ball = MetricBall(1 / unit(eltype(rate)))
  ExponentialTransiogram(rate, ball)
end

"""
    ratematrix(t)

Return the transition rate matrix of the transiogram `t`.
"""
ratematrix(t::ExponentialTransiogram) = t.rate

(t::ExponentialTransiogram)(h) = exp(h * t.rate)
