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

metricball(c::Covariance) = metricball(c.γ)

(c::Covariance)(h) = sill(c.γ) - c.γ(h)

"""
    isstationary(c)

Check if covariance `c` possesses the 2nd-order stationary property.
"""
isstationary(c::Covariance) = isstationary(c.γ)

"""
    sill(c)

Return the sill of the covariance `c`.
"""
sill(c::Covariance) = sill(c.γ)

"""
    nugget(c)

Return the nugget of the covariance `c`.
"""
nugget(c::Covariance) = nugget(c.γ)

"""
    scale(c, s)

Scale metric ball of covariance `c` with strictly
positive scaling factor `s`.
"""
scale(c::C, s::Real) where {C<:Covariance} = C(scale(c.γ, s))

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, c::T) where {T<:Covariance}
  name = string(nameof(T))
  _showcompact(io, name, c.γ)
end

function Base.show(io::IO, ::MIME"text/plain", c::T) where {T<:Covariance}
  name = string(nameof(T))
  _showfull(io, name, c.γ)
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
