# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PiecewiseLinearTransiogram(ball,  abscissas, ordinates)

A piecewise-linear transiogram model with given metric `ball`,
`abscissas` and matrix `ordinates`.

    PiecewiseLinearTransiogram(abscissas, ordinates; range)
    PiecewiseLinearTransiogram(abscissas, ordinates; ranges, rotation)

Alternatively, build metric `ball` from single `range` or from
multiple `ranges` and `rotation` matrix.

## References

* Li, W. & Zhang, C. 2005. [Application of Transiograms to Markov Chain
  Simulation and Spatial Uncertainty Assessment of Land-Cover Classes]
  (https://www.tandfonline.com/doi/abs/10.2747/1548-1603.42.4.297)

* Li, W. & Zhang, C. 2010. [Linear interpolation and joint model fitting
  of experimental transiograms for Markov chain simulation of categorical
  spatial variables](https://www.tandfonline.com/doi/abs/10.1080/13658810903127991)
"""
struct PiecewiseLinearTransiogram{B<:MetricBall,ℒ<:Len,M} <: Transiogram
  ball::B
  abscissas::Vector{ℒ}
  ordinates::Vector{M}

  function PiecewiseLinearTransiogram{B,ℒ,M}(ball, abscissas, ordinates) where {B<:MetricBall,ℒ<:Len,M}
    @assert all(h -> h ≥ zero(h), abscissas) "invalid abscissas for piecewise linear model"
    @assert all(allequal, eachcol(last(ordinates))) "invalid ordinates for piecewise linear model"
    new{B,ℒ,M}(ball, abscissas, ordinates)
  end
end

function PiecewiseLinearTransiogram(ball::MetricBall, abscissas::AbstractVector, ordinates::AbstractVector)
  # metric ball
  b = ball

  # abscissa vector
  x = float(aslen.(abscissas))

  # ordinate matrices
  Y = [SMatrix{size(Yᵢ)...}(Yᵢ) for Yᵢ in ordinates]

  # proportions from last ordinate matrix
  Yₙ = last(Y)
  yₙ = diag(Yₙ)
  k = length(yₙ)
  p = if all(iszero, yₙ)
    k⁻¹ = eltype(yₙ)(1 / k)
    SVector{k}(k⁻¹ for i in 1:k)
  else
    normalize(yₙ, 1)
  end

  # replace last ordinate matrix based on proportions
  Y[end] = SMatrix{k,k}(p[j] for i in 1:k, j in 1:k)

  PiecewiseLinearTransiogram{typeof(b),eltype(x),eltype(Y)}(b, x, Y)
end

PiecewiseLinearTransiogram(
  abscissas::AbstractVector,
  ordinates::AbstractVector;
  range=1.0,
  ranges=nothing,
  rotation=I
) = PiecewiseLinearTransiogram(rangeball(range, ranges, rotation), abscissas, ordinates)

constructor(::PiecewiseLinearTransiogram) = PiecewiseLinearTransiogram

Base.range(t::PiecewiseLinearTransiogram) = maximum(meanlengths(t))

scale(t::PiecewiseLinearTransiogram, s::Real) = PiecewiseLinearTransiogram(s * t.ball, t.abscissas, t.ordinates)

function meanlengths(t::PiecewiseLinearTransiogram)
  h₁ = t.abscissas[1]
  h₂ = t.abscissas[2]
  t₁ = diag(t.ordinates[1])
  t₂ = diag(t.ordinates[2])
  a = (t₂ - t₁) / (h₂ - h₁)
  b = t₁ - h₁ * a
  l = -(b ./ a)
  r = ustrip(maximum(radii(t.ball)))
  l .* r
end

proportions(t::PiecewiseLinearTransiogram) = Tuple(diag(last(t.ordinates)))

function (t::PiecewiseLinearTransiogram)(h)
  h′ = aslen(h)
  hs = t.abscissas
  if h′ < first(hs) # left extrapolation
    ((first(hs) - h′) * I + h′ * first(t.ordinates)) / first(hs)
  elseif h′ ≥ last(hs) # right extrapolation
    last(t.ordinates)
  else # middle interpolation
    k = findlast(≤(h′), hs)
    ((hs[k + 1] - h′) * t.ordinates[k] + (h′ - hs[k]) * t.ordinates[k + 1]) / (hs[k + 1] - hs[k])
  end
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, t::PiecewiseLinearTransiogram)
  ioctx = IOContext(io, :compact => true)
  summary(ioctx, t)
  print(ioctx, "(")
  if hasequalradii(t.ball)
    print(ioctx, "range: ", first(radii(t.ball)))
  else
    print(ioctx, "ranges: ", radii(t.ball))
  end
  print(ioctx, ")")
end

function Base.show(io::IO, ::MIME"text/plain", t::PiecewiseLinearTransiogram)
  ioctx = IOContext(io, :compact => true, :limit => true)
  summary(ioctx, t)
  println(ioctx)
  if hasequalradii(t.ball)
    println(ioctx, "├─ range: ", first(radii(t.ball)))
  else
    println(ioctx, "├─ ranges: ", radii(t.ball))
  end
  print(ioctx, "├─ abscissas: ")
  _printlnvec(ioctx, t.abscissas, 3)
  print(ioctx, "├─ ordinates: ")
  _printvec(ioctx, t.ordinates, 3)
end
