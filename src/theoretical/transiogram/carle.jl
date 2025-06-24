# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CarleTransiogram(ratex, ratey, ratez, ...)

A transiogram with transition rate matrices `ratex`, `ratey`, `ratez`, ...
along axes of a Cartesian coordinate reference system.

## Examples

```julia
# lengths and proportions
l = (1u"m", 2u"m", 3u"m")
p = (1/3, 1/3, 1/3)

# horizontal to vertical anisotropy ratio
ratex = GeoStatsFunctions.baseratematrix(10 *. l, p)
ratey = GeoStatsFunctions.baseratematrix(10 *. l, p)
ratez = GeoStatsFunctions.baseratematrix(1  *. l, p)

# Carle's transiogram model
CarleTransiogram(ratex, ratey, ratez)
```

See also [`MatrixExponentialTransiogram`](@ref).

## References

* Carle, S.F. & Fogg, G.E. 1996. [Transition probability-based
  indicator geostatistics](https://link.springer.com/article/10.1007/BF02083656)

* Carle et al 1998. [Conditional Simulation of Hydrofacies Architecture:
  A Transition Probability/Markov Approach](https://doi.org/10.2110/sepmcheg.01.147)
"""
struct CarleTransiogram{N,R<:StaticMatrix} <: Transiogram
  rates::NTuple{N,R}

  function CarleTransiogram{N,R}(rates) where {N,R<:StaticMatrix}
    @assert allequal(size(rate) for rate in rates) "transition rate matrices must have equal size"
    @assert allequal(size(first(rates))) "transition rate matrices must be square"
    new(rates)
  end
end

function CarleTransiogram(rates::Tuple)
  srates = ntuple(i -> SMatrix{size(rates[i])...}(rates[i]), length(rates))
  urates = ntuple(i -> asinvlen.(srates[i]), length(srates))
  CarleTransiogram{length(urates),eltype(urates)}(urates)
end

CarleTransiogram(rates::AbstractMatrix...) = CarleTransiogram(rates)

constructor(::CarleTransiogram) = CarleTransiogram

function metricball(t::CarleTransiogram)
  R̃ = t.rates
  N = length(R̃)
  j = argmax(meanlengths(t))
  MetricBall(ntuple(i -> 1 / -diag(R̃[i])[j], N))
end

Base.range(t::CarleTransiogram) = maximum(meanlengths(t))

scale(t::CarleTransiogram, s::Real) = CarleTransiogram(ntuple(i -> t.rates[i] / s, length(t.rates)))

function meanlengths(t::CarleTransiogram)
  R̃ = t.rates
  N = length(R̃)
  l = ntuple(N) do i
    1 ./ -diag(R̃[i])
  end
  Tuple(max.(l...))
end

proportions(t::CarleTransiogram) = Tuple(normalize(diag(t(100range(t))), 1))

function (t::CarleTransiogram)(p₁::Point, p₂::Point)
  h⃗ = p₂ - p₁
  h = norm(h⃗)
  R̃ = t.rates
  k = size(first(R̃), 1)

  # handle corner case 
  iszero(h) && return SMatrix{k,k}(one(h) * I(k))

  # rate matrix along direction h⃗
  r(i, j) = √sum(abs2, h⃗[n] * R̃[n][i, j] for n in eachindex(h⃗))
  R = SMatrix{k,k}(r(i, j) for i in 1:k, j in 1:k) / h

  # normalized transition matrix
  T = exp(h * R)
  T ./ sum(T, dims=2)
end

function (t::CarleTransiogram)(h)
  R̃ = t.rates
  N = length(R̃)
  j = argmax(i -> maximum(1 ./ -diag(R̃[i])), 1:N)
  p₁ = Point(ntuple(i -> zero(h), N))
  p₂ = Point(ntuple(i -> i == j ? h : zero(h), N))
  t(p₁, p₂)
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, t::CarleTransiogram)
  ioctx = IOContext(io, :compact => true)
  summary(ioctx, t)
  print(ioctx, "(meanlengths: ", meanlengths(t), ", proportions: ", proportions(t), ")")
end

function Base.show(io::IO, ::MIME"text/plain", t::CarleTransiogram)
  ioctx = IOContext(io, :compact => true, :limit => true)
  summary(ioctx, t)
  println(ioctx)
  print(ioctx, "├─ rates: ")
  _printvec(ioctx, collect(t.rates), 3)
end
