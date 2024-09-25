# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Variogram

A theoretical variogram model (e.g. Gaussian variogram).
"""
abstract type Variogram end

"""
    constructor(γ)

Return the type constructor of the variogram `γ`.
"""
function constructor end

"""
    isstationary(γ)

Check if variogram `γ` possesses the 2nd-order stationary property.
"""
isstationary(γ::Variogram) = isstationary(typeof(γ))

"""
    isisotropic(γ)

Tells whether or not variogram `γ` is isotropic.
"""
isisotropic(γ::Variogram) = isisotropic(γ.ball)

"""
    sill(γ)

Return the sill of the variogram `γ`.
"""
sill(γ::Variogram) = γ.sill

"""
    nugget(γ)

Return the nugget of the variogram `γ`.
"""
nugget(γ::Variogram) = γ.nugget

"""
    metricball(γ)

Return the metric ball of the variogram `γ`.
"""
metricball(γ::Variogram) = γ.ball

"""
    range(γ)

Return the maximum range of the variogram `γ`.
"""
Base.range(γ::Variogram) = maximum(radii(γ.ball))

"""
    scale(γ, s)

Scale metric ball of variogram `γ` with strictly
positive scaling factor `s`.
"""
function scale(γ::Variogram, s::Real)
  V = constructor(γ)
  V(s * metricball(γ); sill=sill(γ), nugget=nugget(γ))
end

"""
    structures(γ)

Return the individual structures of a (possibly nested)
variogram as a tuple. The structures are the total nugget
`cₒ`, and the coefficients (or contributions) `c[i]` for the
remaining non-trivial structures `g[i]` after normalization
(i.e. sill=1, nugget=0).

## Examples

```julia
γ₁ = GaussianVariogram(nugget=1, sill=2)
γ₂ = SphericalVariogram(nugget=2, sill=3)

γ = 2γ₁ + 3γ₂

cₒ, c, g = structures(γ)
```
"""
function structures(γ::Variogram)
  cₒ = nugget(γ)
  c = sill(γ) - nugget(γ)
  T = typeof(c)
  γ = @set γ.sill = one(T)
  γ = @set γ.nugget = zero(T)
  cₒ, (c,), (γ,)
end

"""
    γ(u, v)

Evaluate the variogram at points `u` and `v`.
"""
function (γ::Variogram)(u::Point, v::Point)
  d = metric(γ.ball)
  h = evaluate(d, u, v)
  γ(h)
end

"""
    γ(U, v)

Evaluate the variogram at geometry `U` and point `v`.
"""
function (γ::Variogram)(U::Geometry, v::Point)
  us = variosample(γ, U)
  mean(γ(u, v) for u in us)
end

"""
    γ(u, V)

Evaluate the variogram at point `u` and geometry `V`.
"""
(γ::Variogram)(u::Point, V::Geometry) = γ(V, u)

"""
    γ(U, V)

Evaluate the variogram at geometries `U` and `V`.
"""
function (γ::Variogram)(U::Geometry, V::Geometry)
  us = variosample(γ, U)
  vs = variosample(γ, V)
  mean(γ(u, v) for u in us, v in vs)
end

"""
    pairwise(γ, domain)

Evaluate variogram `γ` between all elements in the `domain`.
"""
function pairwise(γ::Variogram, domain)
  u = first(domain)
  n = length(domain)
  R = returntype(γ, u, u)
  Γ = Matrix{R}(undef, n, n)
  pairwise!(Γ, γ, domain)
end

"""
    pairwise!(Γ, γ, domain)

Evaluates covariance `γ` between all elements in the `domain` in-place, filling the matrix `Γ`.
"""
function pairwise!(Γ, γ::Variogram, domain)
  n = length(domain)
  @inbounds for j in 1:n
    vⱼ = domain[j]
    sⱼ = variosample(γ, vⱼ)
    for i in (j + 1):n
      vᵢ = domain[i]
      sᵢ = variosample(γ, vᵢ)
      Γ[i, j] = mean(γ(pᵢ, pⱼ) for pᵢ in sᵢ, pⱼ in sⱼ)
    end
    Γ[j, j] = mean(γ(pⱼ, pⱼ) for pⱼ in sⱼ, pⱼ in sⱼ)
    for i in 1:(j - 1)
      Γ[i, j] = Γ[j, i] # leverage the symmetry
    end
  end
  Γ
end

"""
    pairwise(γ, domain₁, domain₂)

Evaluate variogram `γ` between all elements of `domain₁` and `domain₂`.
"""
function pairwise(γ::Variogram, domain₁, domain₂)
  u = first(domain₁)
  v = first(domain₂)
  m = length(domain₁)
  n = length(domain₂)
  R = returntype(γ, u, v)
  Γ = Matrix{R}(undef, m, n)
  pairwise!(Γ, γ, domain₁, domain₂)
end

"""
    pairwise!(Γ, γ, domain₁, domain₂)

Evaluates covariance `γ` between all elements of `domain₁` and `domain₂` in-place, filling the matrix `Γ`.
"""
function pairwise!(Γ, γ::Variogram, domain₁, domain₂)
  m = length(domain₁)
  n = length(domain₂)
  @inbounds for j in 1:n
    vⱼ = domain₂[j]
    sⱼ = variosample(γ, vⱼ)
    for i in 1:m
      vᵢ = domain₁[i]
      sᵢ = variosample(γ, vᵢ)
      Γ[i, j] = mean(γ(pᵢ, pⱼ) for pᵢ in sᵢ, pⱼ in sⱼ)
    end
  end
  Γ
end

"""
    returntype(γ, u, v)

Return type of γ(u, v).
"""
returntype(γ::Variogram, u, v) = typeof(γ(u, v))

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, γ::T) where {T<:Variogram}
  name = string(nameof(T))
  _showcompact(io, name, γ)
end

function Base.show(io::IO, ::MIME"text/plain", γ::T) where {T<:Variogram}
  name = string(nameof(T))
  _showfull(io, name, γ)
end

function _showcompact(io, name, γ::T) where {T<:Variogram}
  params = String[]
  for fn in fieldnames(T)
    val = getfield(γ, fn)
    if val isa MetricBall
      if isisotropic(val)
        r = first(radii(val))
        push!(params, "range: $r")
      else
        r = Tuple(radii(val))
        push!(params, "ranges: $r")
      end
      d = nameof(typeof(metric(val)))
      push!(params, "distance: $d")
    else
      push!(params, "$fn: $val")
    end
  end
  print(io, name, "(", join(params, ", "), ")")
end

function _showfull(io, name, γ::T) where {T<:Variogram}
  header = isisotropic(γ) ? name : name * " (anisotropic)"
  params = String[]
  fnames = fieldnames(T)
  len = length(fnames)
  for (i, fn) in enumerate(fnames)
    div = i == len ? "└─" : "├─"
    val = getfield(γ, fn)
    if val isa MetricBall
      if isisotropic(val)
        r = first(radii(val))
        push!(params, "├─ range: $r")
      else
        r = Tuple(radii(val))
        push!(params, "└─ ranges: $r")
      end
      m = nameof(typeof(metric(val)))
      push!(params, "$div distance: $m")
    else
      push!(params, "$div $(fn): $val")
    end
  end
  println(io, header)
  print(io, join(params, "\n"))
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("variogram/gaussian.jl")
include("variogram/exponential.jl")
include("variogram/spherical.jl")
include("variogram/matern.jl")
include("variogram/cubic.jl")
include("variogram/pentaspherical.jl")
include("variogram/sinehole.jl")
include("variogram/power.jl")
include("variogram/nugget.jl")
include("variogram/circular.jl")
