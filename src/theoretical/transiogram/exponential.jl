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

"""
    ratematrix(t)

Return the transition rate matrix of the transiogram `t`.
"""
ratematrix(t::ExponentialTransiogram) = t.rate

levels(t::ExponentialTransiogram) = t.levs

(t::ExponentialTransiogram)(h) = exp(h * t.rate)
