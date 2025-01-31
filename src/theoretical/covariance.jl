# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Covariance

A theoretical covariance function (e.g. Gaussian covariance).
"""
abstract type Covariance <: GeoStatsFunction end

# ---------------------
# GEOSTATSFUNCTION API
# ---------------------

isstationary(cov::Covariance) = isstationary(cov.γ)

isbanded(::Type{<:Covariance}) = true

metricball(cov::Covariance) = metricball(cov.γ)

(cov::Covariance)(h) = sill(cov.γ) - cov.γ(h)

# ---------------
# COVARIANCE API
# ---------------

"""
    sill(cov)

Return the sill of the covariance `cov`.
"""
sill(cov::Covariance) = sill(cov.γ)

"""
    nugget(cov)

Return the nugget of the covariance `cov`.
"""
nugget(cov::Covariance) = nugget(cov.γ)

"""
    structures(cov)

Return the individual structures of a (possibly composite)
covariance as a tuple. The structures are the total nugget,
and the coefficients (or contributions) for for the remaining
non-trivial structures after normalization (i.e. sill=1, nugget=0).

## Examples

```julia
cov₁ = GaussianCovariance(nugget=1, sill=2)
cov₂ = SphericalCovariance(nugget=2, sill=3)

structures(2cov₁ + 3cov₂)
```
"""
function structures(cov::C) where {C<:Covariance}
  cₒ, c, γ = structures(cov.γ)
  cₒ, c, C.(γ)
end

"""
    scale(cov, s)

Scale metric ball of covariance `cov` with strictly
positive scaling factor `s`.
"""
scale(cov::C, s::Real) where {C<:Covariance} = C(scale(cov.γ, s))

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, cov::Covariance)
  ioctx = IOContext(io, :compact => true)
  summary(ioctx, cov)
  print(ioctx, "(")
  _printfields(ioctx, cov.γ, singleline=true)
  print(ioctx, ")")
end

function Base.show(io::IO, ::MIME"text/plain", cov::Covariance)
  ioctx = IOContext(io, :compact => true, :limit => true)
  summary(ioctx, cov)
  println(ioctx)
  _printfields(ioctx, cov.γ)
end

# heper macro to define covariances
macro defcov(CovType, VarioType)
  docstring = """
      $(CovType)(args..., kwargs...)

  A covariance function derived from the corresponding variogram function.

  Please see [`$(VarioType)`](@ref) for available parameters.
  """
  expr = quote
    @doc $docstring struct $CovType{V<:$VarioType} <: Covariance
      γ::V
    end

    $CovType(args...; kwargs...) = $CovType($VarioType(args...; kwargs...))
  end
  esc(expr)
end

@defcov CircularCovariance CircularVariogram

@defcov CubicCovariance CubicVariogram

@defcov ExponentialCovariance ExponentialVariogram

@defcov GaussianCovariance GaussianVariogram

@defcov MaternCovariance MaternVariogram

@defcov PentaSphericalCovariance PentaSphericalVariogram

@defcov SineHoleCovariance SineHoleVariogram

@defcov SphericalCovariance SphericalVariogram
