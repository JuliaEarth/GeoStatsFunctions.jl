# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

include("fitting/algorithms.jl")

# -----------------
# GEOSTATSFUNCTION
# -----------------

"""
    fit(F, f, algo=WeightedLeastSquares(); kwargs...)

Fit theoretical geostatistical function of type `F` to empirical function `f` using algorithm `algo`.

Optionally fix theoretical parameters like `range`, `sill` and `nugget` in the `kwargs`.

## Examples

```julia
julia> fit(SphericalVariogram, g)
julia> fit(ExponentialVariogram, g)
julia> fit(ExponentialVariogram, g, sill=1.0)
julia> fit(ExponentialVariogram, g, maxsill=1.0)
julia> fit(GaussianVariogram, g, WeightedLeastSquares())
```
"""
fit(F::Type{<:GeoStatsFunction}, f::EmpiricalGeoStatsFunction, algo::FitAlgo=WeightedLeastSquares(); kwargs...) =
  _fit(F, f, algo; kwargs...) |> first

"""
    fit(Fs, f, algo=WeightedLeastSquares(); kwargs...)

Fit theoretical geostatistical functions of types `Fs` to empirical function `f`
using algorithm `algo` and return the one with minimum error.

## Examples

```julia
julia> fit([SphericalVariogram, ExponentialVariogram], g)
```
"""
function fit(Fs, f::EmpiricalGeoStatsFunction, algo::FitAlgo=WeightedLeastSquares(); kwargs...)
  # fit each variogram type
  res = [_fit(F, f, algo; kwargs...) for F in Fs]
  fs, ϵs = first.(res), last.(res)

  # return best candidate
  fs[argmin(ϵs)]
end

"""
    fit(F, f, weightfun; kwargs...)

Convenience method that forwards the weighting function `weightfun`
to the `WeightedLeastSquares` algorithm.

## Examples

```julia
fit(SphericalVariogram, g, h -> exp(-h))
fit(Variogram, g, h -> exp(-h/100))
```
"""
fit(F, f::EmpiricalGeoStatsFunction, weightfun::Function; kwargs...) =
  fit(F, f, WeightedLeastSquares(weightfun); kwargs...)

# ----------
# VARIOGRAM
# ----------

"""
    fit(Variogram, g, algo=WeightedLeastSquares(); kwargs...)

Fit all stationary `Variogram` models to empirical variogram `g`
using algorithm `algo` and return the one with minimum error.

## Examples

```julia
julia> fit(Variogram, g)
julia> fit(Variogram, g, h -> 1 / h)
```
"""
function fit(::Type{Variogram}, g::EmpiricalVariogram, algo::FitAlgo=WeightedLeastSquares(); kwargs...)
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
  fit(Gs, g, algo; kwargs...)
end

"""
    fit(Transiogram, t, algo=WeightedLeastSquares(); kwargs...)

Fit all theoretical `Transiogram` models to empirical transiogram `t`
using algorithm `algo` and return the one with minimum error.

## Examples

```julia
julia> fit(Transiogram, t)
julia> fit(Transiogram, t, h -> 1 / h)
```
"""
function fit(::Type{Transiogram}, t::EmpiricalTransiogram, algo::FitAlgo=WeightedLeastSquares(); kwargs...)
  Ts = (
    ExponentialTransiogram,
    GaussianTransiogram,
    LinearTransiogram,
    MatrixExponentialTransiogram,
    SphericalTransiogram
  )
  fit(Ts, t, algo; kwargs...)
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

function _optimize(J, L, λ, l, u, θₒ)
  sol = Optim.optimize(θ -> J(θ) + λ * L(θ), l, u, θₒ)
  ϵ = Optim.minimum(sol)
  θ = Optim.minimizer(sol)
  θ, ϵ
end
