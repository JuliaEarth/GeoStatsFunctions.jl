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
struct CarleTransiogram{N,R<:StaticMatrix,P<:NTuple} <: Transiogram
  rates::NTuple{N,R}
  proportions::P

  function CarleTransiogram{N,R,P}(rates, proportions) where {N,R<:StaticMatrix,P<:NTuple}
    @assert allequal(size(rate) for rate in rates) "transition rate matrices must have equal size"
    @assert allequal(size(first(rates))) "transition rate matrices must be square"
    @assert all(p -> 0 ≤ p ≤ 1, proportions) "proportions must be in [0, 1] interval"
    @assert sum(proportions) ≈ 1 "proportions must sum up to one"
    new(rates, proportions)
  end
end

function CarleTransiogram(rates::Tuple)
  # convert rates to unitful static matrices
  srates = ntuple(i -> SMatrix{size(rates[i])...}(rates[i]), length(rates))
  urates = ntuple(i -> asinvlen.(srates[i]), length(srates))

  # proportions from first rate matrix
  R = first(urates)
  r = maximum(1 ./ -diag(R))
  props = Tuple(diag(exp(100r * R)))

  CarleTransiogram{length(urates),eltype(urates),typeof(props)}(urates, props)
end

CarleTransiogram(rates::AbstractMatrix...) = CarleTransiogram(rates)

function CarleTransiogram()
  R = baseratematrix((1.0, 1.0), (0.5, 0.5))
  CarleTransiogram(R, R, R)
end

constructor(::CarleTransiogram) = CarleTransiogram

function isisotropic(t::CarleTransiogram)
  R = t.rates
  N = length(R)
  l = ntuple(N) do i
    1 ./ -diag(R[i])
  end
  allequal(l)
end

metricball(t::CarleTransiogram) = MetricBall(1u"m")

Base.range(t::CarleTransiogram) = maximum(meanlengths(t))

scale(t::CarleTransiogram, s::Real) = CarleTransiogram(ntuple(i -> t.rates[i] / s, length(t.rates)))

function meanlengths(t::CarleTransiogram)
  R = t.rates
  N = length(R)
  l = ntuple(N) do i
    1 ./ -diag(R[i])
  end
  Tuple(max.(l...))
end

proportions(t::CarleTransiogram) = t.proportions

function (t::CarleTransiogram)(p₁::Point, p₂::Point)
  # lag vector and norm
  h⃗ = p₂ - p₁
  h = norm(h⃗)

  # number of levels and proportions
  R = t.rates
  k = size(first(R), 1)
  p = t.proportions

  # handle corner case 
  iszero(h) && return SMatrix{k,k}(one(h) * I(k))

  # Eq. 22 of Carle et al 1998
  w(i, j, d) = h⃗[d] ≥ zero(h⃗[d]) ? R[d][i, j] : (p[j] / p[i]) * R[d][j, i]
  r(i, j) = ((i == j) ? -1 : 1) * √sum(abs2, h⃗[d] * w(i, j, d) for d in eachindex(h⃗))

  # Eq. 20 and 21 of Carle et al 1998
  A = SMatrix{k - 1,k - 1}(r(i, j) for i in 1:(k - 1), j in 1:(k - 1)) / h
  b = -sum(p[1:(k - 1)] .* A, dims=1) / p[k]
  c = -sum([A; b], dims=2)
  Rₕ = [[A; b] c]

  # transition matrix
  exp(h * Rₕ)
end

function (t::CarleTransiogram)(h)
  R = t.rates
  N = length(R)
  j = argmax(i -> maximum(1 ./ -diag(R[i])), 1:N)
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
