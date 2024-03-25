# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Covariance

Parent type of all covariance functions (e.g. Gaussian covariance).
"""
abstract type Covariance end

"""
    cov(x₁, x₂)

Evaluate the covariance at objects `x₁` and `x₁`.
"""
(cov::Covariance)(x₁, x₂) = sill(cov.γ) - cov.γ(x₁, x₂)

"""
    isstationary(cov)

Check if covariance `cov` possesses the 2nd-order stationary property.
"""
isstationary(cov::Covariance) = isstationary(cov.γ)

"""
    isisotropic(cov)

Tells whether or not covariance `cov` is isotropic.
"""
isisotropic(cov::Covariance) = isisotropic(cov.γ)

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
    metricball(cov)

Return the metric ball of the covariance `cov`.
"""
metricball(cov::Covariance) = metricball(cov.γ)

"""
    range(cov)

Return the maximum range of the covariance `cov`."""
Base.range(cov::Covariance) = range(cov.γ)

"""
    scale(cov, s)

Scale metric ball of covariance `cov` with strictly
positive scaling factor `s`.
"""
scale(cov::Cov, s::Real) where {Cov <: Covariance} = Cov(scale(cov.γ, s))

"""
    pairwise(cov, domain)

Evaluate covariance `cov` between all elements in the `domain`.

    pairwise(cov, domain₁, domain₂)

Evaluate covariance `cov` between all elements of `domain₁` and `domain₂`.
"""
pairwise(cov::Covariance, args...) = sill(cov.γ) .- pairwise(cov.γ, args...)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, cov::T) where {T<:Covariance}
  name = string(nameof(T))
  _showcompact(io, name, cov.γ)
end

function Base.show(io::IO, ::MIME"text/plain", cov::T) where {T<:Covariance}
  name = string(nameof(T))
  _showfull(io, name, cov.γ)
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

@defcov PentasphericalCovariance PentasphericalVariogram

@defcov SineHoleCovariance SineHoleVariogram

@defcov SphericalCovariance SphericalVariogram
