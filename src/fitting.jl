# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    fit(F, f, [w]; parameters..., maxparameters...)

Fit theoretical geostatistical function of type `F` to empirical function `f`.
The weighting function `w` maps lags to importance weights. It is optional and
defaults to the bin counts stored in each lag of the empirical function `f`.

Theoretical `parameters` like `range`, `sill` and `nugget` as well as their
maximum values like `maxrange`, `maxsill` and `maxnugget` can be fixed as
keyword arguments. They are converted into optimization constraints in the
weighted least squares fitting algorithm.

The type `F` can be abstract like `Variogram` or `Transiogram`, in which case
all fittable (e.g., stationary) subtypes are fitted and the one with minimum
error is returned. See [`fitany`](@ref) for custom lists of types.

## Examples

```julia
julia> fit(SphericalVariogram, g)
julia> fit(ExponentialVariogram, g, range=1.0)
julia> fit(GaussianVariogram, g, sill=1.0, nugget=0.1)
julia> fit(ExponentialVariogram, g, maxsill=1.0)
julia> fit(SphericalVariogram, g, h -> 1/h)
julia> fit(Variogram, g, range=1.0)
julia> fit(Transiogram, t, h -> 1/h)
```
"""
fit(F::Type{<:GeoStatsFunction}, f::EmpiricalGeoStatsFunction, w=nothing; kwargs...) = _fit(F, f, w; kwargs...) |> first

function fit(::Type{Variogram}, g::EmpiricalVariogram, w=nothing; kwargs...)
  Gs = (
    CircularVariogram,
    CubicVariogram,
    ExponentialVariogram,
    GaussianVariogram,
    MaternVariogram,
    PentaSphericalVariogram,
    SineHoleVariogram,
    SphericalVariogram
  )
  fitany(Gs, g, w; kwargs...)
end

function fit(::Type{Transiogram}, t::EmpiricalTransiogram, w=nothing; kwargs...)
  Ts =
    (ExponentialTransiogram, GaussianTransiogram, LinearTransiogram, MatrixExponentialTransiogram, SphericalTransiogram)
  fitany(Ts, t, w; kwargs...)
end

"""
    fitany([F₁, F₂, ...], f, [w]; parameters..., maxparameters...)

Fit theoretical geostatistical functions of types `F₁, F₂, ...`
to empirical function `f` and return the one with minimum error.

## Examples

```julia
julia> fitany([SphericalVariogram, ExponentialVariogram], g)
```

See [`fit`](@ref) for details on the arguments and keyword arguments.
"""
function fitany(Fs, f::EmpiricalGeoStatsFunction, w=nothing; kwargs...)
  sols = [_fit(F, f, w; kwargs...) for F in Fs]
  f, _ = argmin(last, sols)
  f
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("fitting/variograms.jl")
include("fitting/transiograms.jl")

# -----------------
# HELPER FUNCTIONS
# -----------------

_ustrip(u, x) = x
_ustrip(u, x::Quantity) = ustrip(u, x)

_weights(f, x, n) = isnothing(f) ? n / sum(n) : map(xᵢ -> ustrip(f(xᵢ)), x)

function _optimize(J, θₗ, θᵤ, θₒ)
  s = Optim.optimize(J, θₗ, θᵤ, θₒ, LBFGSB())
  ϵ = Optim.minimum(s)
  θ = Optim.minimizer(s)
  θ, ϵ
end
