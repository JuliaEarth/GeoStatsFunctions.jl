# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExponentialTransiogram(R)

An exponential transiogram with transition rate matrix `R`.

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)
"""
struct ExponentialTransiogram{M<:StaticMatrix} <: Transiogram
  R::M

  function ExponentialTransiogram{M}(R) where {M<:StaticMatrix}
    if !allequal(size(R))
      throw(ArgumentError("transition rate matrix must be square"))
    end
    new(R)
  end
end

function ExponentialTransiogram(R::AbstractMatrix)
  S = SMatrix{size(R)...}(R)
  ExponentialTransiogram{typeof(S)}(S)
end

"""
    ratematrix(t)

Return the transition rate matrix of the transiogram `t`.
"""
ratematrix(t::ExponentialTransiogram) = t.R
