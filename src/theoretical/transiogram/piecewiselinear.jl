# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PiecewiseLinearTransiogram(abscissas, ordinates)

A piecewise-linear transiogram model with `abscissas` and
matrix `ordinates` obtained from an [`EmpiricalTransiogram`](@ref).

    PiecewiseLinearTransiogram(ball, abscissas, ordinates)

Alternatively, use a custom metric `ball`.

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
  ordinfinity::M
end

function PiecewiseLinearTransiogram(ball::MetricBall, abscissas::AbstractVector, ordinates::AbstractMatrix)
  n = length(abscissas)
  m = size(ordinates, 1)

  # metric ball
  b = ball

  # abscissa vector
  x = float(aslen.(abscissas))

  # ordinate matrices
  Y = map(1:n) do k
    SMatrix{m,m}(ordinates[i, j][k] for i in 1:m, j in 1:m)
  end

  # ordinate matrix at infinity
  yₖ = diag(last(Y))
  p = if all(iszero, yₖ)
    # corner case with uniform proportions
    m⁻¹ = eltype(yₖ)(1 / m)
    SVector{m}(m⁻¹ for i in 1:m)
  else
    normalize(yₖ, 1)
  end
  Y∞ = SMatrix{m,m}(p[j] for i in 1:m, j in 1:m)

  PiecewiseLinearTransiogram(b, x, Y, Y∞)
end

function PiecewiseLinearTransiogram(abscissas::AbstractVector, ordinates::AbstractMatrix)
  ball = MetricBall(oneunit(eltype(abscissas)))
  PiecewiseLinearTransiogram(ball, abscissas, ordinates)
end

proportions(t::PiecewiseLinearTransiogram) = Tuple(diag(t.ordinfinity))

function (t::PiecewiseLinearTransiogram)(h)
  h′ = aslen(h)
  hs = t.abscissas
  if h′ < first(hs) # left extrapolation
    ((first(hs) - h′) * I + h′ * first(t.ordinates)) / first(hs)
  elseif h′ > last(hs) # right extrapolation
    t.ordinfinity
  else # middle interpolation
    k = findlast(<(h′), hs)
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
